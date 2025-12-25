#ifndef SRC_NMIE_BATCH_HPP_
#define SRC_NMIE_BATCH_HPP_

#include <vector>
#include <complex>
#include <algorithm>
#include "nmie-precision.hpp"
#include "nmie-basic.hpp"

namespace nmie {

struct MieBatchInput {
  std::vector<double> x;
  std::vector<std::complex<double>> m;
};

struct MieBatchOutput {
  std::vector<double> Qext;
  std::vector<double> Qsca;
  std::vector<double> Qabs;
  std::vector<double> Qbk;
  std::vector<double> Qpr;
  std::vector<double> g;
  std::vector<double> Albedo;
};

template <typename FloatType, typename Engine = HighwayEngine<FloatType>>
MieBatchOutput RunMieBatch(const MieBatchInput& input) {
  MieBatchOutput output;
  size_t N = input.x.size();
  if (input.m.size() != N) throw std::runtime_error("Input size mismatch");

  output.Qext.resize(N);
  output.Qsca.resize(N);
  output.Qabs.resize(N);
  output.Qbk.resize(N);
  output.Qpr.resize(N);
  output.g.resize(N);
  output.Albedo.resize(N);

  const hn::ScalableTag<FloatType> d;
  const size_t lanes = hn::Lanes(d);

  for (size_t i = 0; i < N; i += lanes) {
    // 1. Load Batch
    std::vector<FloatType> x_batch(lanes, 1.0); // Avoid x=0 division
    std::vector<FloatType> m_re_batch(lanes, 1.0);
    std::vector<FloatType> m_im_batch(lanes, 0.0);
    std::vector<FloatType> nmax_vals(lanes);
    
    size_t remaining = N - i;
    size_t current_batch_size = std::min(lanes, remaining);

    FloatType max_x = 0;
    int max_nmax = 0;

    for (size_t j = 0; j < current_batch_size; ++j) {
      x_batch[j] = input.x[i + j];
      m_re_batch[j] = input.m[i + j].real();
      m_im_batch[j] = input.m[i + j].imag();
      
      FloatType x = x_batch[j];
      if (x > max_x) max_x = x;
      
      int nmax = std::round(x + 11 * std::pow(x, (1.0 / 3.0)) + 1);
      nmax_vals[j] = static_cast<FloatType>(nmax);
      if (nmax > max_nmax) max_nmax = nmax;
    }
    // Pad the rest with valid but dummy values (e.g. x=1, m=1) to avoid NaNs/Inf in calculation
    // The results for these padded lanes will be ignored.
    for (size_t j = current_batch_size; j < lanes; ++j) {
        nmax_vals[j] = 0; // Should stop immediately
    }

    auto vx = hn::Load(d, x_batch.data());
    auto vmr = hn::Load(d, m_re_batch.data());
    auto vmi = hn::Load(d, m_im_batch.data());
    auto vnmax = hn::Load(d, nmax_vals.data());

    typename Engine::ComplexV mL_v = { vmr, vmi };
    typename Engine::ComplexV z_v = { hn::Mul(vx, vmr), hn::Mul(vx, vmi) };
    typename Engine::ComplexV x_real_v = { vx, hn::Zero(d) };

    // 2. Calculate Special Functions for the whole batch
    int calc_nmax = max_nmax + 15; // Buffer for stability
    std::vector<typename Engine::ComplexV> D1_z(calc_nmax + 1), D1_x(calc_nmax + 1), 
                                           Psi_x(calc_nmax + 1), Zeta_x(calc_nmax + 1), D3_x(calc_nmax + 1);
    std::vector<typename Engine::ComplexV> PsiZeta_x(calc_nmax + 1);

    evalDownwardD1<FloatType, Engine>(z_v, D1_z);
    evalDownwardD1<FloatType, Engine>(x_real_v, D1_x);
    evalUpwardPsi<FloatType, Engine>(x_real_v, D1_x, Psi_x);
    evalUpwardD3<FloatType, Engine>(x_real_v, D1_x, D3_x, PsiZeta_x);

    // Efficiency factors accumulation registers
    auto vQext = hn::Zero(d);
    auto vQsca = hn::Zero(d);
    auto vQbktmp_re = hn::Zero(d);
    auto vQbktmp_im = hn::Zero(d);
    auto vQpr_sum = hn::Zero(d);

    typename Engine::ComplexV an_prev = { hn::Zero(d), hn::Zero(d) };
    typename Engine::ComplexV bn_prev = { hn::Zero(d), hn::Zero(d) };
    
    // Convergence threshold
    auto vThreshold = hn::Set(d, 1e-14);

    // 3. Sum loop (Horizontal SIMD: across particles)
    for (int n = 1; n < calc_nmax; ++n) {
      auto vn = hn::Set(d, static_cast<FloatType>(n));
      auto mask = hn::Le(vn, vnmax);
      
      // Optimization: Check if all active lanes have finished (n > nmax)
      // But nmax varies. We can check if mask is all false.
      // However, we also want to check convergence of an/bn.
      
      auto Zeta_n = Engine::div(PsiZeta_x[n], Psi_x[n]);
      auto Zeta_nm1 = Engine::div(PsiZeta_x[n-1], Psi_x[n-1]);

      auto an = calc_an<FloatType, Engine>(n, vx, D1_z[n], mL_v, 
                                        Psi_x[n], Zeta_n, Psi_x[n-1], Zeta_nm1);
      auto bn = calc_bn<FloatType, Engine>(n, vx, D1_z[n], mL_v, 
                                        Psi_x[n], Zeta_n, Psi_x[n-1], Zeta_nm1);

      // Mask out contributions beyond nmax for each lane
      an.re = hn::IfThenElse(mask, an.re, hn::Zero(d));
      an.im = hn::IfThenElse(mask, an.im, hn::Zero(d));
      bn.re = hn::IfThenElse(mask, bn.re, hn::Zero(d));
      bn.im = hn::IfThenElse(mask, bn.im, hn::Zero(d));
      
      // Dynamic Convergence Check
      auto an_mag = hn::Add(hn::Abs(an.re), hn::Abs(an.im)); // Approx magnitude
      auto bn_mag = hn::Add(hn::Abs(bn.re), hn::Abs(bn.im));
      auto max_mag = hn::Max(an_mag, bn_mag);
      
      auto converged = hn::Lt(max_mag, vThreshold);
      auto finished = hn::Not(mask);
      auto done = hn::Or(converged, finished);
      
      if (hn::AllTrue(d, done)) break;

      if (n > max_nmax) break;

      // mult = (2n + 1)
      auto n_float = vn;
      auto mult = hn::Add(hn::Add(n_float, n_float), hn::Set(d, 1.0));
      
      // Qext += (2n+1) * Re(an + bn)
      vQext = hn::Add(vQext, hn::Mul(mult, hn::Add(an.re, bn.re)));

      // Qsca += (2n+1) * (|an|^2 + |bn|^2)
      auto an_mag2 = hn::Add(hn::Mul(an.re, an.re), hn::Mul(an.im, an.im));
      auto bn_mag2 = hn::Add(hn::Mul(bn.re, bn.re), hn::Mul(bn.im, bn.im));
      vQsca = hn::Add(vQsca, hn::Mul(mult, hn::Add(an_mag2, bn_mag2)));

      // Qbk sum: Qbktmp += (2n+1) * (-1)^n * (an - bn)
      // Scalar code uses (1.0 - 2.0 * (n1 % 2)) where n1=n.
      // n=1 -> -1. n=2 -> +1.
      // My previous code: n=1 -> +1.
      // So I flip the signs to match Scalar.
      auto sign_n = ((n + 1) % 2 == 0) ? hn::Set(d, -1.0) : hn::Set(d, 1.0);
      vQbktmp_re = hn::Add(vQbktmp_re, hn::Mul(hn::Mul(mult, sign_n), hn::Sub(an.re, bn.re)));
      vQbktmp_im = hn::Add(vQbktmp_im, hn::Mul(hn::Mul(mult, sign_n), hn::Sub(an.im, bn.im)));

      // Qpr interaction term (n-1 and n)
      if (n > 1) {
        auto nm1 = hn::Sub(n_float, hn::Set(d, 1.0));
        auto pre1 = hn::Div(hn::Mul(nm1, hn::Add(n_float, hn::Set(d, 1.0))), n_float);
        // Re(a_{n-1}*a_n* + b_{n-1}*b_n*)
        auto re_aa = hn::Add(hn::Mul(an_prev.re, an.re), hn::Mul(an_prev.im, an.im));
        auto re_bb = hn::Add(hn::Mul(bn_prev.re, bn.re), hn::Mul(bn_prev.im, bn.im));
        vQpr_sum = hn::Add(vQpr_sum, hn::Mul(pre1, hn::Add(re_aa, re_bb)));
      }
      // Qpr internal term
      auto pre2 = hn::Div(mult, hn::Mul(n_float, hn::Add(n_float, hn::Set(d, 1.0))));
      auto re_ab = hn::Add(hn::Mul(an.re, bn.re), hn::Mul(an.im, bn.im));
      vQpr_sum = hn::Add(vQpr_sum, hn::Mul(pre2, re_ab));

      an_prev = an; bn_prev = bn;
    }

    // Final normalization: 2/x^2
    auto x2 = hn::Mul(vx, vx);
    auto norm = hn::Div(hn::Set(d, 2.0), x2);
    vQext = hn::Mul(vQext, norm);
    vQsca = hn::Mul(vQsca, norm);
    
    // Qbk = |Qbktmp|^2 / x^2
    auto qbk_mag2 = hn::Add(hn::Mul(vQbktmp_re, vQbktmp_re), hn::Mul(vQbktmp_im, vQbktmp_im));
    auto vQbk = hn::Div(qbk_mag2, x2);

    // Qpr = Qext - 4/x^2 * sum
    auto norm_pr = hn::Div(hn::Set(d, 4.0), x2);
    auto vQpr = hn::Sub(vQext, hn::Mul(norm_pr, vQpr_sum));

    // g = (Qext - Qpr) / Qsca
    auto vg = hn::Div(hn::Sub(vQext, vQpr), vQsca);
    // Handle Qsca=0 case? If Qsca=0, g=0.
    vg = hn::IfThenElse(hn::Eq(vQsca, hn::Zero(d)), hn::Zero(d), vg);

    // Albedo = Qsca / Qext
    auto vAlbedo = hn::Div(vQsca, vQext);
    vAlbedo = hn::IfThenElse(hn::Eq(vQext, hn::Zero(d)), hn::Zero(d), vAlbedo);

    // Store results
    std::vector<FloatType> res_ext(lanes), res_sca(lanes), res_bk(lanes), res_pr(lanes), res_g(lanes), res_alb(lanes);
    hn::Store(vQext, d, res_ext.data());
    hn::Store(vQsca, d, res_sca.data());
    hn::Store(vQbk, d, res_bk.data());
    hn::Store(vQpr, d, res_pr.data());
    hn::Store(vg, d, res_g.data());
    hn::Store(vAlbedo, d, res_alb.data());

    for (size_t j = 0; j < current_batch_size; ++j) {
      output.Qext[i + j] = res_ext[j];
      output.Qsca[i + j] = res_sca[j];
      output.Qabs[i + j] = res_ext[j] - res_sca[j];
      output.Qbk[i + j] = res_bk[j];
      output.Qpr[i + j] = res_pr[j];
      output.g[i + j] = res_g[j];
      output.Albedo[i + j] = res_alb[j];
    }
  }
  return output;
}

} // namespace nmie

#endif // SRC_NMIE_BATCH_HPP_
