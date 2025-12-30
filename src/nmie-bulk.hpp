#ifndef SRC_NMIE_BULK_HPP_
#define SRC_NMIE_BULK_HPP_

#include <vector>
#include <complex>
#include <algorithm>
#include "nmie.hpp"
#include "special-functions-impl.hpp"
#include "nmie-basic.hpp"

namespace nmie {

template <typename FloatType>
class BulkSphere {
public:
    BulkSphere(const std::vector<std::vector<FloatType>>& x,
               const std::vector<std::vector<std::complex<FloatType>>>& m,
               int nmax = -1) 
        : x_(x), m_(m), nmax_(nmax) {
        if (x.empty() || m.empty()) throw std::invalid_argument("Empty input");
        L_ = x.size();
        num_spheres_ = x[0].size();
        // Validate dimensions
        for (const auto& layer : x) if (layer.size() != num_spheres_) throw std::invalid_argument("Mismatch in x dimensions");
        if (m.size() != L_) throw std::invalid_argument("Mismatch in m layers");
        for (const auto& layer : m) if (layer.size() != num_spheres_) throw std::invalid_argument("Mismatch in m dimensions");
        
        Qext.resize(num_spheres_);
        Qsca.resize(num_spheres_);
        Qabs.resize(num_spheres_);
        Qbk.resize(num_spheres_);
        Qpr.resize(num_spheres_);
        g.resize(num_spheres_);
        Albedo.resize(num_spheres_);
    }

    void RunMieCalculation() {
        using Engine = HighwayEngine<FloatType>;
        const size_t lanes = Engine::Lanes();
        
        MieBuffers<FloatType, Engine> buffers;
        std::vector<std::complex<FloatType>> an_vec, bn_vec;

        for (size_t i = 0; i < num_spheres_; i += lanes) {
            size_t current_batch_size = std::min(lanes, num_spheres_ - i);
            
            // Determine nmax for this batch
            int batch_nmax = nmax_;
            if (batch_nmax == -1) {
                FloatType max_x = 0;
                for (size_t j = 0; j < current_batch_size; ++j) {
                    if (x_[L_-1][i+j] > max_x) max_x = x_[L_-1][i+j];
                }
                batch_nmax = static_cast<int>(std::round(max_x + 11 * std::pow(max_x, 1.0/3.0) + 16));
            }

            buffers.resize(batch_nmax, L_);
            an_vec.resize((batch_nmax + 1) * lanes);
            bn_vec.resize((batch_nmax + 1) * lanes);
            
            auto load_val = [&](const std::vector<FloatType>& v, size_t offset) {
                alignas(64) FloatType tmp[64] = {0};
                for(size_t k=0; k<current_batch_size; ++k) tmp[k] = v[offset+k];
                if (current_batch_size < lanes) {
                     for(size_t k=current_batch_size; k<lanes; ++k) tmp[k] = 1.0; 
                }
                return Engine::load(tmp);
            };
            
            auto load_complex = [&](const std::vector<std::complex<FloatType>>& v, size_t offset) {
                alignas(64) FloatType re[64] = {0};
                alignas(64) FloatType im[64] = {0};
                for(size_t k=0; k<current_batch_size; ++k) {
                    re[k] = v[offset+k].real();
                    im[k] = v[offset+k].imag();
                }
                if (current_batch_size < lanes) {
                     for(size_t k=current_batch_size; k<lanes; ++k) { re[k] = 1.0; im[k] = 0.0; }
                }
                return Engine::make_complex(Engine::load(re), Engine::load(im));
            };

            auto get_x = [&](int l) { return load_val(x_[l], i); };
            auto get_m = [&](int l) { return load_complex(m_[l], i); };

            calcScattCoeffsKernel<FloatType, Engine>(
                batch_nmax, L_, -1, get_x, get_m,
                buffers, an_vec, bn_vec
            );
            
            // Calculate vnmax
            alignas(64) FloatType nmax_arr[64];
            for (size_t j = 0; j < current_batch_size; ++j) {
                FloatType x_val = x_[L_-1][i+j];
                nmax_arr[j] = std::round(x_val + 11 * std::pow(x_val, 1.0/3.0) + 16);
            }
            if (current_batch_size < lanes) {
                for(size_t k=current_batch_size; k<lanes; ++k) nmax_arr[k] = 0.0;
            }
            auto vnmax = Engine::load(nmax_arr);
            
            auto sum_qext = Engine::set(0.0);
            auto sum_qsca = Engine::set(0.0);
            auto sum_qbk_re = Engine::set(0.0);
            auto sum_qbk_im = Engine::set(0.0);
            
            for (int n = 0; n < batch_nmax; ++n) {
                auto n_val = Engine::set(static_cast<FloatType>(n+1));
                auto active_mask = Engine::le(n_val, vnmax);
                
                auto an_val = Engine::load(&an_vec[n*lanes]);
                auto bn_val = Engine::load(&bn_vec[n*lanes]);
                
                // Sanitize NaNs
                auto an_re_raw = Engine::get_real(an_val);
                auto an_im_raw = Engine::get_imag(an_val);
                auto bn_re_raw = Engine::get_real(bn_val);
                auto bn_im_raw = Engine::get_imag(bn_val);
                
                auto nan_mask_an = Engine::lor(Engine::neq(an_re_raw, an_re_raw), Engine::neq(an_im_raw, an_im_raw));
                auto nan_mask_bn = Engine::lor(Engine::neq(bn_re_raw, bn_re_raw), Engine::neq(bn_im_raw, bn_im_raw));
                
                auto zero_c = Engine::make_complex(Engine::set(0.0), Engine::set(0.0));
                an_val = Engine::select(nan_mask_an, zero_c, an_val);
                bn_val = Engine::select(nan_mask_bn, zero_c, bn_val);
                
                // Store back sanitized values
                Engine::store(an_val, &an_vec[n*lanes]);
                Engine::store(bn_val, &bn_vec[n*lanes]);
                
                an_val = Engine::select(active_mask, an_val, zero_c);
                bn_val = Engine::select(active_mask, bn_val, zero_c);
                
                auto two_n_plus_1 = Engine::set(static_cast<FloatType>(2*(n+1) + 1));
                
                auto an_re = Engine::get_real(an_val);
                auto an_im = Engine::get_imag(an_val);
                auto bn_re = Engine::get_real(bn_val);
                auto bn_im = Engine::get_imag(bn_val);
                
                // Qext += (2n+1) * Re(an + bn)
                sum_qext = Engine::add(sum_qext, Engine::mul(two_n_plus_1, Engine::add(an_re, bn_re)));
                
                // Qsca += (2n+1) * (|an|^2 + |bn|^2)
                auto an_sq = Engine::add(Engine::mul(an_re, an_re), Engine::mul(an_im, an_im));
                auto bn_sq = Engine::add(Engine::mul(bn_re, bn_re), Engine::mul(bn_im, bn_im));
                sum_qsca = Engine::add(sum_qsca, Engine::mul(two_n_plus_1, Engine::add(an_sq, bn_sq)));
                
                // Qbk += (2n+1) * (-1)^n * (an - bn)
                FloatType sign_val = ((n+1) % 2 == 1) ? -1.0 : 1.0;
                auto sign = Engine::set(sign_val);
                
                auto diff_re = Engine::sub(an_re, bn_re);
                auto diff_im = Engine::sub(an_im, bn_im);
                
                sum_qbk_re = Engine::add(sum_qbk_re, Engine::mul(sign, Engine::mul(two_n_plus_1, diff_re)));
                sum_qbk_im = Engine::add(sum_qbk_im, Engine::mul(sign, Engine::mul(two_n_plus_1, diff_im)));
            }
            
            auto xL = load_val(x_[L_-1], i);
            auto xL_sq = Engine::mul(xL, xL);
            auto factor = Engine::div(Engine::set(2.0), xL_sq);
            
            auto qext = Engine::mul(factor, sum_qext);
            auto qsca = Engine::mul(factor, sum_qsca);
            auto qabs = Engine::sub(qext, qsca);
            
            auto qbk_sq = Engine::add(Engine::mul(sum_qbk_re, sum_qbk_re), Engine::mul(sum_qbk_im, sum_qbk_im));
            auto qbk = Engine::mul(Engine::div(factor, Engine::set(2.0)), qbk_sq);
            
            // Store results
            alignas(64) FloatType qext_arr[64];
            alignas(64) FloatType qsca_arr[64];
            alignas(64) FloatType qabs_arr[64];
            alignas(64) FloatType qbk_arr[64];
            
            Engine::store(qext, qext_arr);
            Engine::store(qsca, qsca_arr);
            Engine::store(qabs, qabs_arr);
            Engine::store(qbk, qbk_arr);
            
            for(size_t k=0; k<current_batch_size; ++k) {
                Qext[i+k] = qext_arr[k];
                Qsca[i+k] = qsca_arr[k];
                Qabs[i+k] = qabs_arr[k];
                Qbk[i+k] = qbk_arr[k];
            }
        }
    }

    // Getters
    std::vector<FloatType> Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;

private:
    std::vector<std::vector<FloatType>> x_; // [layer][sphere]
    std::vector<std::vector<std::complex<FloatType>>> m_; // [layer][sphere]
    int nmax_;
    size_t L_;
    size_t num_spheres_;
};

}
#endif
