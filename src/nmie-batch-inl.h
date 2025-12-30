namespace hwy {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;
using nmie::MieBatchInput;
using nmie::MieBatchOutput;
using nmie::evalDownwardD1;
using nmie::evalUpwardPsi;
using nmie::evalUpwardD3;
using nmie::calc_an;
using nmie::calc_bn;

// Define the SIMD Engine inside the namespace so it uses the correct target-specific types
#include "highway-engine.hpp"

template <typename FloatType>
MieBatchOutput RunMieBatchImpl(const MieBatchInput& input) {
    using Engine = HighwayEngine<FloatType>;
    MieBatchOutput output;
    size_t N = input.x.size();
    
    output.Qext.resize(N); output.Qsca.resize(N); output.Qabs.resize(N);
    output.Qbk.resize(N); output.Qpr.resize(N); output.g.resize(N);
    output.Albedo.resize(N);

    size_t num_angles = input.theta.size();
    if (num_angles > 0) {
      output.S1.resize(N);
      output.S2.resize(N);
      for (size_t i = 0; i < N; ++i) {
        output.S1[i].resize(num_angles);
        output.S2[i].resize(num_angles);
      }
    }

    const hn::ScalableTag<FloatType> d;
    const size_t lanes = hn::Lanes(d);

    // Calculate global max nmax
    FloatType global_max_x = 0;
    for (auto val : input.x) {
        if (val > global_max_x) global_max_x = val;
    }
    int global_nmax = static_cast<int>(std::round(global_max_x + 11 * std::pow(global_max_x, 1.0/3.0) + 16));

    nmie::MieBuffers<FloatType, Engine> buffers;
    buffers.resize(global_nmax, 1);
    
    std::vector<std::complex<FloatType>> an_vec, bn_vec;
    an_vec.reserve((global_nmax + 1) * lanes);
    bn_vec.reserve((global_nmax + 1) * lanes);

    for (size_t i = 0; i < N; i += lanes) {
        size_t current_batch_size = std::min(lanes, N - i);
        alignas(64) FloatType xb[64] = {0}, mrb[64] = {0}, mib[64] = {0}, nb[64] = {0};
        
        FloatType max_x = 0;
        for (size_t j = 0; j < current_batch_size; ++j) {
            xb[j] = input.x[i + j];
            mrb[j] = input.m[i + j].real();
            mib[j] = input.m[i + j].imag();
            if (xb[j] > max_x) max_x = xb[j];
            nb[j] = static_cast<FloatType>(std::round(xb[j] + 11 * std::pow(xb[j], (1.0 / 3.0)) + 1));
        }
        
        // SIMD loop logic (Horizontal across particles)
        auto vx = hn::LoadU(d, xb);
        auto vmr = hn::LoadU(d, mrb);
        auto vmi = hn::LoadU(d, mib);
        auto vnmax = hn::LoadU(d, nb);

        typename Engine::ComplexV mL_v = { vmr, vmi };

        int calc_nmax = static_cast<int>(std::round(max_x + 11 * std::pow(max_x, 1.0/3.0) + 16));
        
        an_vec.resize((calc_nmax + 1) * lanes);
        bn_vec.resize((calc_nmax + 1) * lanes);

        // Getters for the kernel
        // For batch processing of single spheres, we treat it as L=1.
        // get_x(0) returns the x vector for the batch.
        // get_m(0) returns the m vector for the batch.
        auto get_x = [&](int /*l*/) { return vx; };
        auto get_m = [&](int /*l*/) { return mL_v; };

        nmie::calcScattCoeffsKernel<FloatType, Engine>(
            calc_nmax, 1, -1, get_x, get_m,
            buffers, an_vec, bn_vec
        );

        auto vQext = hn::Zero(d), vQsca = hn::Zero(d), vQpr_sum = hn::Zero(d);
        typename Engine::ComplexV vQbk = {hn::Zero(d), hn::Zero(d)};
        
        std::vector<typename Engine::ComplexV> vS1(num_angles, {hn::Zero(d), hn::Zero(d)});
        std::vector<typename Engine::ComplexV> vS2(num_angles, {hn::Zero(d), hn::Zero(d)});
        
        nmie::sumMieSeriesKernel<FloatType, Engine>(
            calc_nmax,
            vnmax,
            an_vec.data(),
            bn_vec.data(),
            input.theta,
            vQext,
            vQsca,
            vQpr_sum,
            vQbk,
            vS1,
            vS2
        );

        auto x2 = hn::Mul(vx, vx);
        auto norm = hn::Div(hn::Set(d, 2.0), x2);
        
        alignas(64) FloatType r_ext[64], r_sca[64], r_bk[64], r_pr[64];
        hn::Store(hn::Mul(vQext, norm), d, r_ext);
        hn::Store(hn::Mul(vQsca, norm), d, r_sca);
        
        auto vQbk_mag_sq = hn::Add(hn::Mul(vQbk.re, vQbk.re), hn::Mul(vQbk.im, vQbk.im));
        hn::Store(hn::Div(vQbk_mag_sq, x2), d, r_bk);
        
        hn::Store(hn::Sub(hn::Mul(vQext, norm), hn::Mul(hn::Div(hn::Set(d, 4.0), x2), vQpr_sum)), d, r_pr);

        for (size_t j = 0; j < current_batch_size; ++j) {
            output.Qext[i + j] = r_ext[j];
            output.Qsca[i + j] = r_sca[j];
            output.Qabs[i + j] = r_ext[j] - r_sca[j];
            output.Qbk[i + j] = r_bk[j];
            output.Qpr[i + j] = r_pr[j];
            output.g[i + j] = (r_sca[j] > 1e-12) ? (r_ext[j] - r_pr[j]) / r_sca[j] : 0.0;
            output.Albedo[i + j] = (r_ext[j] > 1e-12) ? r_sca[j] / r_ext[j] : 0.0;

            if (num_angles > 0) {
              for (size_t k = 0; k < num_angles; ++k) {
                alignas(64) FloatType r_S1_re[64], r_S1_im[64], r_S2_re[64],
                    r_S2_im[64];
                hn::Store(vS1[k].re, d, r_S1_re);
                hn::Store(vS1[k].im, d, r_S1_im);
                hn::Store(vS2[k].re, d, r_S2_re);
                hn::Store(vS2[k].im, d, r_S2_im);

                output.S1[i + j][k] =
                    std::complex<double>(r_S1_re[j], r_S1_im[j]);
                output.S2[i + j][k] =
                    std::complex<double>(r_S2_re[j], r_S2_im[j]);
              }
            }
        }
    }
    return output;
}

// Concrete double implementation for exporting
MieBatchOutput RunMieBatchDouble(const MieBatchInput& input) {
    return RunMieBatchImpl<double>(input);
}

}  // namespace HWY_NAMESPACE
}  // namespace hwy
