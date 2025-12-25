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

    const hn::ScalableTag<FloatType> d;
    const size_t lanes = hn::Lanes(d);

    for (size_t i = 0; i < N; i += lanes) {
        size_t current_batch_size = std::min(lanes, N - i);
        alignas(64) FloatType xb[64], mrb[64], mib[64], nb[64];
        
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
        typename Engine::ComplexV z_v = { hn::Mul(vx, vmr), hn::Mul(vx, vmi) };
        typename Engine::ComplexV x_real_v = { vx, hn::Zero(d) };

        int calc_nmax = static_cast<int>(std::round(max_x + 11 * std::pow(max_x, 1.0/3.0) + 16));
        
        std::vector<typename Engine::ComplexV> D1_z(calc_nmax + 1), D1_x(calc_nmax + 1), D3_x(calc_nmax + 1),
                                                Psi_x(calc_nmax + 1), PsiZeta_x(calc_nmax + 1);

        evalDownwardD1<FloatType, Engine>(z_v, D1_z);
        evalDownwardD1<FloatType, Engine>(x_real_v, D1_x);
        evalUpwardPsi<FloatType, Engine>(x_real_v, D1_x, Psi_x);
        evalUpwardD3<FloatType, Engine>(x_real_v, D1_x, D3_x, PsiZeta_x);

        auto vQext = hn::Zero(d), vQsca = hn::Zero(d), vQbktmp_re = hn::Zero(d), vQbktmp_im = hn::Zero(d), vQpr_sum = hn::Zero(d);
        typename Engine::ComplexV an_prev = { hn::Zero(d), hn::Zero(d) }, bn_prev = { hn::Zero(d), hn::Zero(d) };

        for (int n = 1; n < calc_nmax; ++n) {
            auto vn = hn::Set(d, static_cast<FloatType>(n));
            auto mask = hn::Le(vn, vnmax);
            if (hn::AllTrue(d, hn::Not(mask))) break;

            auto Zeta_n = Engine::div(PsiZeta_x[n], Psi_x[n]);
            auto Zeta_nm1 = Engine::div(PsiZeta_x[n-1], Psi_x[n-1]);

            auto an = calc_an<FloatType, Engine>(n, vx, D1_z[n], mL_v, Psi_x[n], Zeta_n, Psi_x[n-1], Zeta_nm1);
            auto bn = calc_bn<FloatType, Engine>(n, vx, D1_z[n], mL_v, Psi_x[n], Zeta_n, Psi_x[n-1], Zeta_nm1);

            auto mult = hn::Add(hn::Add(vn, vn), hn::Set(d, 1.0));
            vQext = hn::Add(vQext, hn::IfThenElse(mask, hn::Mul(mult, hn::Add(an.re, bn.re)), hn::Zero(d)));
            vQsca = hn::Add(vQsca, hn::IfThenElse(mask, hn::Mul(mult, hn::Add(hn::Mul(an.re, an.re), hn::Add(hn::Mul(an.im, an.im), hn::Add(hn::Mul(bn.re, bn.re), hn::Mul(bn.im, bn.im))))), hn::Zero(d)));
            
            auto sign_n = ((n + 1) % 2 == 0) ? hn::Set(d, -1.0) : hn::Set(d, 1.0);
            vQbktmp_re = hn::Add(vQbktmp_re, hn::IfThenElse(mask, hn::Mul(hn::Mul(mult, sign_n), hn::Sub(an.re, bn.re)), hn::Zero(d)));
            vQbktmp_im = hn::Add(vQbktmp_im, hn::IfThenElse(mask, hn::Mul(hn::Mul(mult, sign_n), hn::Sub(an.im, bn.im)), hn::Zero(d)));

            if (n > 1) {
                auto re_aa = hn::Add(hn::Mul(an_prev.re, an.re), hn::Mul(an_prev.im, an.im));
                auto re_bb = hn::Add(hn::Mul(bn_prev.re, bn.re), hn::Mul(bn_prev.im, bn.im));
                vQpr_sum = hn::Add(vQpr_sum, hn::IfThenElse(mask, hn::Mul(hn::Div(hn::Mul(hn::Sub(vn, hn::Set(d, 1.0)), hn::Add(vn, hn::Set(d, 1.0))), vn), hn::Add(re_aa, re_bb)), hn::Zero(d)));
            }
            vQpr_sum = hn::Add(vQpr_sum, hn::IfThenElse(mask, hn::Mul(hn::Div(mult, hn::Mul(vn, hn::Add(vn, hn::Set(d, 1.0)))), hn::Add(hn::Mul(an.re, bn.re), hn::Mul(an.im, bn.im))), hn::Zero(d)));

            an_prev = an; bn_prev = bn;
        }

        auto x2 = hn::Mul(vx, vx);
        auto norm = hn::Div(hn::Set(d, 2.0), x2);
        
        alignas(64) FloatType r_ext[64], r_sca[64], r_bk[64], r_pr[64];
        hn::Store(hn::Mul(vQext, norm), d, r_ext);
        hn::Store(hn::Mul(vQsca, norm), d, r_sca);
        hn::Store(hn::Div(hn::Add(hn::Mul(vQbktmp_re, vQbktmp_re), hn::Mul(vQbktmp_im, vQbktmp_im)), x2), d, r_bk);
        hn::Store(hn::Sub(hn::Mul(vQext, norm), hn::Mul(hn::Div(hn::Set(d, 4.0), x2), vQpr_sum)), d, r_pr);

        for (size_t j = 0; j < current_batch_size; ++j) {
            output.Qext[i + j] = r_ext[j];
            output.Qsca[i + j] = r_sca[j];
            output.Qabs[i + j] = r_ext[j] - r_sca[j];
            output.Qbk[i + j] = r_bk[j];
            output.Qpr[i + j] = r_pr[j];
            output.g[i + j] = (r_sca[j] > 1e-12) ? (r_ext[j] - r_pr[j]) / r_sca[j] : 0.0;
            output.Albedo[i + j] = (r_ext[j] > 1e-12) ? r_sca[j] / r_ext[j] : 0.0;
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
