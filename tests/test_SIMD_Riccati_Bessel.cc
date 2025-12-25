#include "gtest/gtest.h"
#include "hwy/highway.h"
#include "../src/special-functions-impl.hpp"
#include "../src/nmie-basic.hpp"
#include "../src/nmie-precision.hpp"
#include "test_spec_functions_data.hpp"

// #define HWY_TARGET_INCLUDE "tests/test_SIMD_Riccati_Bessel.cc"
// #include "hwy/foreach_target.h"

namespace nmie {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

TEST(SIMDRiccatiBessel, ComplexCotMatchScalar) {
  const hn::ScalableTag<double> d;
  using Engine = HighwayEngine<double>;
  
  // Test values (from Yang case or similar)
  double test_re = 1.05 * 80.0;
  double test_im = 1.0 * 80.0;
  
  // 1. Scalar Reference
  std::complex<double> z_scalar(test_re, test_im);
  auto ref_scalar = complex_cot<double, ScalarEngine>(z_scalar);

  // 2. SIMD Execution
  typename Engine::ComplexV z_simd{hn::Set(d, test_re), hn::Set(d, test_im)};
  auto res_simd = complex_cot<double, Engine>(z_simd);
  
  // 3. Verification
  std::vector<double> out_re(hn::Lanes(d));
  std::vector<double> out_im(hn::Lanes(d));
  hn::Store(res_simd.re, d, out_re.data());
  hn::Store(res_simd.im, d, out_im.data());

  for (size_t i = 0; i < hn::Lanes(d); ++i) {
    EXPECT_NEAR(out_re[i], ref_scalar.real(), 1e-13);
    EXPECT_NEAR(out_im[i], ref_scalar.imag(), 1e-13);
  }
}

TEST(SIMDRiccatiBessel, D1RecurrenceMatchScalar) {
  const hn::ScalableTag<double> d;
  using Engine = HighwayEngine<double>;
  
  double test_re = 1.05 * 80.0;
  double test_im = 1.0 * 80.0;
  int nmax = 10; // Small nmax for test
  
  // 1. Scalar Reference
  std::complex<double> z_scalar(test_re, test_im);
  std::vector<std::complex<double>> D1_scalar(nmax + 1);
  evalDownwardD1<double, ScalarEngine>(z_scalar, D1_scalar);

  // 2. SIMD Execution
  typename Engine::ComplexV z_simd{hn::Set(d, test_re), hn::Set(d, test_im)};
  std::vector<typename Engine::ComplexV> D1_simd(nmax + 1);
  evalDownwardD1<double, Engine>(z_simd, D1_simd);
  
  // 3. Verification
  std::vector<double> out_re(hn::Lanes(d));
  std::vector<double> out_im(hn::Lanes(d));
  
  for (int n = 0; n <= nmax; ++n) {
    hn::Store(D1_simd[n].re, d, out_re.data());
    hn::Store(D1_simd[n].im, d, out_im.data());
    
    for (size_t i = 0; i < hn::Lanes(d); ++i) {
      EXPECT_NEAR(out_re[i], D1_scalar[n].real(), 1e-13) << "Mismatch at n=" << n;
      EXPECT_NEAR(out_im[i], D1_scalar[n].imag(), 1e-13) << "Mismatch at n=" << n;
    }
  }
}

TEST(SIMDRiccatiBessel, D3RecurrenceMatchScalar) {
  const hn::ScalableTag<double> d;
  using Engine = HighwayEngine<double>;
  
  double test_re = 1.05 * 80.0;
  double test_im = 1.0 * 80.0;
  int nmax = 10;
  
  // 1. Scalar Reference
  std::complex<double> z_scalar(test_re, test_im);
  std::vector<std::complex<double>> D1_scalar(nmax + 1);
  std::vector<std::complex<double>> D3_scalar(nmax + 1);
  std::vector<std::complex<double>> PsiZeta_scalar(nmax + 1);
  
  evalDownwardD1<double, ScalarEngine>(z_scalar, D1_scalar);
  evalUpwardD3<double, ScalarEngine>(z_scalar, D1_scalar, D3_scalar, PsiZeta_scalar);

  // 2. SIMD Execution
  typename Engine::ComplexV z_simd{hn::Set(d, test_re), hn::Set(d, test_im)};
  std::vector<typename Engine::ComplexV> D1_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> D3_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> PsiZeta_simd(nmax + 1);
  
  evalDownwardD1<double, Engine>(z_simd, D1_simd);
  evalUpwardD3<double, Engine>(z_simd, D1_simd, D3_simd, PsiZeta_simd);
  
  // 3. Verification
  std::vector<double> out_re(hn::Lanes(d));
  std::vector<double> out_im(hn::Lanes(d));
  
  for (int n = 0; n <= nmax; ++n) {
    // Verify D3
    hn::Store(D3_simd[n].re, d, out_re.data());
    hn::Store(D3_simd[n].im, d, out_im.data());
    for (size_t i = 0; i < hn::Lanes(d); ++i) {
      EXPECT_NEAR(out_re[i], D3_scalar[n].real(), 1e-13) << "D3 Mismatch at n=" << n;
      EXPECT_NEAR(out_im[i], D3_scalar[n].imag(), 1e-13) << "D3 Mismatch at n=" << n;
    }
    
    // Verify PsiZeta
    hn::Store(PsiZeta_simd[n].re, d, out_re.data());
    hn::Store(PsiZeta_simd[n].im, d, out_im.data());
    for (size_t i = 0; i < hn::Lanes(d); ++i) {
      EXPECT_NEAR(out_re[i], PsiZeta_scalar[n].real(), 1e-13) << "PsiZeta Mismatch at n=" << n;
      EXPECT_NEAR(out_im[i], PsiZeta_scalar[n].imag(), 1e-13) << "PsiZeta Mismatch at n=" << n;
    }
  }
}

TEST(SIMDRiccatiBessel, PsiRecurrenceMatchScalar) {
  const hn::ScalableTag<double> d;
  using Engine = HighwayEngine<double>;
  
  double test_re = 1.05 * 80.0;
  double test_im = 1.0 * 80.0;
  int nmax = 10;
  
  // 1. Scalar Reference
  std::complex<double> z_scalar(test_re, test_im);
  std::vector<std::complex<double>> D1_scalar(nmax + 1);
  std::vector<std::complex<double>> Psi_scalar(nmax + 1);
  
  evalDownwardD1<double, ScalarEngine>(z_scalar, D1_scalar);
  evalUpwardPsi<double, ScalarEngine>(z_scalar, D1_scalar, Psi_scalar);

  // 2. SIMD Execution
  typename Engine::ComplexV z_simd{hn::Set(d, test_re), hn::Set(d, test_im)};
  std::vector<typename Engine::ComplexV> D1_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> Psi_simd(nmax + 1);
  
  evalDownwardD1<double, Engine>(z_simd, D1_simd);
  evalUpwardPsi<double, Engine>(z_simd, D1_simd, Psi_simd);
  
  // 3. Verification
  std::vector<double> out_re(hn::Lanes(d));
  std::vector<double> out_im(hn::Lanes(d));
  
  for (int n = 0; n <= nmax; ++n) {
    hn::Store(Psi_simd[n].re, d, out_re.data());
    hn::Store(Psi_simd[n].im, d, out_im.data());
    for (size_t i = 0; i < hn::Lanes(d); ++i) {
      EXPECT_NEAR(out_re[i], Psi_scalar[n].real(), 1e-13) << "Psi Mismatch at n=" << n;
      EXPECT_NEAR(out_im[i], Psi_scalar[n].imag(), 1e-13) << "Psi Mismatch at n=" << n;
    }
  }
}

TEST(SIMDRiccatiBessel, ZetaRecurrenceMatchScalar) {
  const hn::ScalableTag<double> d;
  using Engine = HighwayEngine<double>;
  
  double test_re = 1.05 * 80.0;
  double test_im = 1.0 * 80.0;
  int nmax = 10;
  
  // 1. Scalar Reference
  std::complex<double> z_scalar(test_re, test_im);
  std::vector<std::complex<double>> D1_scalar(nmax + 1);
  std::vector<std::complex<double>> D3_scalar(nmax + 1);
  std::vector<std::complex<double>> PsiZeta_scalar(nmax + 1); // Used for D3 calc
  std::vector<std::complex<double>> Zeta_scalar(nmax + 1);
  
  evalDownwardD1<double, ScalarEngine>(z_scalar, D1_scalar);
  evalUpwardD3<double, ScalarEngine>(z_scalar, D1_scalar, D3_scalar, PsiZeta_scalar);
  evalUpwardZeta<double, ScalarEngine>(z_scalar, D3_scalar, Zeta_scalar);

  // 2. SIMD Execution
  typename Engine::ComplexV z_simd{hn::Set(d, test_re), hn::Set(d, test_im)};
  std::vector<typename Engine::ComplexV> D1_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> D3_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> PsiZeta_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> Zeta_simd(nmax + 1);
  
  evalDownwardD1<double, Engine>(z_simd, D1_simd);
  evalUpwardD3<double, Engine>(z_simd, D1_simd, D3_simd, PsiZeta_simd);
  evalUpwardZeta<double, Engine>(z_simd, D3_simd, Zeta_simd);
  
  // 3. Verification
  std::vector<double> out_re(hn::Lanes(d));
  std::vector<double> out_im(hn::Lanes(d));
  
  for (int n = 0; n <= nmax; ++n) {
    hn::Store(Zeta_simd[n].re, d, out_re.data());
    hn::Store(Zeta_simd[n].im, d, out_im.data());
    for (size_t i = 0; i < hn::Lanes(d); ++i) {
      EXPECT_NEAR(out_re[i], Zeta_scalar[n].real(), 1e-13) << "Zeta Mismatch at n=" << n;
      EXPECT_NEAR(out_im[i], Zeta_scalar[n].imag(), 1e-13) << "Zeta Mismatch at n=" << n;
    }
  }
}

TEST(SIMDRiccatiBessel, D1FullLanesMatchDataset) {
  const hn::ScalableTag<double> d;
  using Engine = HighwayEngine<double>;
  const size_t lanes = hn::Lanes(d);
  
  // We need enough data to fill the lanes
  // We'll take the first 'lanes' entries from an_test_30digits 
  // (which contains x, m, n, result)
  ASSERT_GE(an_test_30digits.size(), lanes);

  std::vector<double> x_vals(lanes), mr_vals(lanes), mi_vals(lanes);
  std::vector<int> n_vals(lanes);
  std::vector<std::complex<double>> expected_scalar(lanes);

  int target_n = 1; // Testing index n=1 for all lanes initially

  for (size_t i = 0; i < lanes; ++i) {
    // an_test_30digits schema: {x, {mr, mi}, n, {res_re, res_im}, err_re, err_im}
    auto& entry = an_test_30digits[i]; 
    x_vals[i] = std::get<0>(entry);
    mr_vals[i] = std::get<1>(entry).real();
    mi_vals[i] = std::get<1>(entry).imag();
    // For D1, we calculate the whole vector up to some nmax
  }

  // Pack into SIMD
  auto vx = hn::Load(d, x_vals.data());
  auto vmr = hn::Load(d, mr_vals.data());
  auto vmi = hn::Load(d, mi_vals.data());

  typename Engine::ComplexV z_simd = { hn::Mul(vx, vmr), hn::Mul(vx, vmi) };

  // Run vectorized D1
  int nmax = 100; 
  std::vector<typename Engine::ComplexV> D1_simd(nmax + 1);
  evalDownwardD1<double, Engine>(z_simd, D1_simd);

  // Verification: Run scalar for each lane and compare
  for (size_t i = 0; i < lanes; ++i) {
    std::complex<double> z_lane(x_vals[i] * mr_vals[i], x_vals[i] * mi_vals[i]);
    std::vector<std::complex<double>> D1_ref(nmax + 1);
    evalDownwardD1<double, ScalarEngine>(z_lane, D1_ref);

    std::vector<double> res_re(lanes), res_im(lanes);
    hn::Store(D1_simd[target_n].re, d, res_re.data());
    hn::Store(D1_simd[target_n].im, d, res_im.data());

    EXPECT_NEAR(res_re[i], D1_ref[target_n].real(), 1e-12) << "Lane " << i << " mismatch at D1[" << target_n << "]";
    EXPECT_NEAR(res_im[i], D1_ref[target_n].imag(), 1e-12) << "Lane " << i << " mismatch at D1[" << target_n << "]";
  }
}

TEST(SIMDRiccatiBessel, AnFullLanesMatchDataset) {
  const hn::ScalableTag<double> d;
  using Engine = HighwayEngine<double>;
  const size_t lanes = hn::Lanes(d);
  
  ASSERT_GE(an_test_30digits.size(), lanes);

  std::vector<double> x_vals(lanes), mr_vals(lanes), mi_vals(lanes);
  std::vector<int> n_vals(lanes);
  
  // We will test for a specific n, say n=1, across different parameters if possible,
  // or just use the same parameters but different lanes to verify consistency.
  // However, the dataset has different n for same x,m.
  // Let's pick entries where n=1 if possible, or just use whatever is there and calculate for max n.
  
  // For simplicity, let's just take the first 'lanes' entries.
  // We will calculate up to nmax=100 (sufficient for small x in dataset)
  int nmax = 100;

  for (size_t i = 0; i < lanes; ++i) {
    auto& entry = an_test_30digits[i]; 
    x_vals[i] = std::get<0>(entry);
    mr_vals[i] = std::get<1>(entry).real();
    mi_vals[i] = std::get<1>(entry).imag();
    n_vals[i] = std::get<2>(entry);
  }

  auto vx = hn::Load(d, x_vals.data());
  auto vmr = hn::Load(d, mr_vals.data());
  auto vmi = hn::Load(d, mi_vals.data());

  // z = x * m
  typename Engine::ComplexV z_simd = { hn::Mul(vx, vmr), hn::Mul(vx, vmi) };
  // XL = x
  typename Engine::RealV XL_simd = vx;
  // mL = m
  typename Engine::ComplexV mL_simd = { vmr, vmi };
  // Ha = 0 (for bulk sphere, Ha=0 usually? No, Ha is related to magnetic permeability or similar boundary conditions?)
  // In standard Mie, for non-magnetic sphere, mu=1.
  // Wait, calc_an signature: calc_an(n, XL, Ha, mL, PsiXL, ZetaXL, PsiXLM1, ZetaXLM1)
  // In nmie-basic.hpp:
  // an_[n] = calc_an(n + 1, x[L - 1], Ha[L - 1][n], m[L - 1], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
  // For a single sphere (bulk), Ha is 0?
  // Let's check calcScattCoeffs in nmie.hpp or nmie-basic.hpp.
  // In nmie-basic.hpp, calcScattCoeffs calls calc_an.
  // It seems Ha is calculated.
  
  // For the purpose of this unit test, we want to verify that calc_an_v produces the same result as calc_an scalar
  // given the SAME inputs. We don't necessarily need to reproduce the full physics of the sphere here,
  // just the mathematical function calc_an.
  
  // So we will generate inputs for calc_an and feed them to both scalar and SIMD versions.
  
  typename Engine::ComplexV Ha_simd;
  
  // We need Psi and Zeta for XL (real).
  // And Ha for z = x*m (complex).
  
  // 1. Calculate Psi and Zeta for XL (real)
  typename Engine::ComplexV z_real_simd = { vx, hn::Zero(d) };
  std::vector<typename Engine::ComplexV> D1_real_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> D3_real_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> PsiZeta_real_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> Psi_real_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> Zeta_real_simd(nmax + 1);

  evalDownwardD1<double, Engine>(z_real_simd, D1_real_simd);
  evalUpwardPsi<double, Engine>(z_real_simd, D1_real_simd, Psi_real_simd);
  evalUpwardD3<double, Engine>(z_real_simd, D1_real_simd, D3_real_simd, PsiZeta_real_simd);
  
  for(int k=0; k<=nmax; ++k) {
      Zeta_real_simd[k] = Engine::div(PsiZeta_real_simd[k], Psi_real_simd[k]);
  }

  // 2. Calculate Ha for z = x*m (complex)
  // Ha corresponds to D1(z) for the core.
  std::vector<typename Engine::ComplexV> D1_complex_simd(nmax + 1);
  evalDownwardD1<double, Engine>(z_simd, D1_complex_simd);
  
  // For single sphere, Ha = D1(z)
  // We need Ha at target_n
  int target_n = 1;
  Ha_simd = D1_complex_simd[target_n];

  // Now calculate an for a specific n, say n=1.
  auto an_simd = calc_an<double, Engine>(target_n, XL_simd, Ha_simd, mL_simd, 
                                         Psi_real_simd[target_n], Zeta_real_simd[target_n], 
                                         Psi_real_simd[target_n-1], Zeta_real_simd[target_n-1]);

  // Verification
  std::vector<double> res_re(lanes), res_im(lanes);
  hn::Store(an_simd.re, d, res_re.data());
  hn::Store(an_simd.im, d, res_im.data());

  for (size_t i = 0; i < lanes; ++i) {
    // Scalar calculation
    std::complex<double> z_lane(x_vals[i] * mr_vals[i], x_vals[i] * mi_vals[i]);
    std::complex<double> z_real_lane(x_vals[i], 0.0);
    std::complex<double> m_lane(mr_vals[i], mi_vals[i]);
    double XL_lane = x_vals[i];
    
    std::vector<std::complex<double>> D1_real_ref(nmax + 1);
    std::vector<std::complex<double>> D3_real_ref(nmax + 1);
    std::vector<std::complex<double>> PsiZeta_real_ref(nmax + 1);
    std::vector<std::complex<double>> Psi_real_ref(nmax + 1);
    std::vector<std::complex<double>> Zeta_real_ref(nmax + 1);

    evalDownwardD1<double, ScalarEngine>(z_real_lane, D1_real_ref);
    evalUpwardPsi<double, ScalarEngine>(z_real_lane, D1_real_ref, Psi_real_ref);
    evalUpwardD3<double, ScalarEngine>(z_real_lane, D1_real_ref, D3_real_ref, PsiZeta_real_ref);
    for(int k=0; k<=nmax; ++k) Zeta_real_ref[k] = PsiZeta_real_ref[k] / Psi_real_ref[k];

    std::vector<std::complex<double>> D1_complex_ref(nmax + 1);
    evalDownwardD1<double, ScalarEngine>(z_lane, D1_complex_ref);
    std::complex<double> Ha_lane = D1_complex_ref[target_n];

    auto an_ref = calc_an<double, ScalarEngine>(target_n, XL_lane, Ha_lane, m_lane,
                                                Psi_real_ref[target_n], Zeta_real_ref[target_n],
                                                Psi_real_ref[target_n-1], Zeta_real_ref[target_n-1]);

    EXPECT_NEAR(res_re[i], an_ref.real(), 1e-12) << "Lane " << i << " mismatch at an[" << target_n << "]";
    EXPECT_NEAR(res_im[i], an_ref.imag(), 1e-12) << "Lane " << i << " mismatch at an[" << target_n << "]";
  }
}

TEST(SIMDRiccatiBessel, BnFullLanesMatchDataset) {
  const hn::ScalableTag<double> d;
  using Engine = HighwayEngine<double>;
  const size_t lanes = hn::Lanes(d);
  
  ASSERT_GE(an_test_30digits.size(), lanes);

  std::vector<double> x_vals(lanes), mr_vals(lanes), mi_vals(lanes);
  int nmax = 100;

  for (size_t i = 0; i < lanes; ++i) {
    auto& entry = an_test_30digits[i]; 
    x_vals[i] = std::get<0>(entry);
    mr_vals[i] = std::get<1>(entry).real();
    mi_vals[i] = std::get<1>(entry).imag();
  }

  auto vx = hn::Load(d, x_vals.data());
  auto vmr = hn::Load(d, mr_vals.data());
  auto vmi = hn::Load(d, mi_vals.data());

  typename Engine::ComplexV z_simd = { hn::Mul(vx, vmr), hn::Mul(vx, vmi) };
  typename Engine::RealV XL_simd = vx;
  typename Engine::ComplexV mL_simd = { vmr, vmi };
  typename Engine::ComplexV Hb_simd; 
  
  // 1. Calculate Psi and Zeta for XL (real)
  typename Engine::ComplexV z_real_simd = { vx, hn::Zero(d) };
  std::vector<typename Engine::ComplexV> D1_real_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> D3_real_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> PsiZeta_real_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> Psi_real_simd(nmax + 1);
  std::vector<typename Engine::ComplexV> Zeta_real_simd(nmax + 1);

  evalDownwardD1<double, Engine>(z_real_simd, D1_real_simd);
  evalUpwardPsi<double, Engine>(z_real_simd, D1_real_simd, Psi_real_simd);
  evalUpwardD3<double, Engine>(z_real_simd, D1_real_simd, D3_real_simd, PsiZeta_real_simd);
  
  for(int k=0; k<=nmax; ++k) {
      Zeta_real_simd[k] = Engine::div(PsiZeta_real_simd[k], Psi_real_simd[k]);
  }

  // 2. Calculate Hb for z = x*m (complex)
  // Hb corresponds to D1(z) for the core (same as Ha for core).
  std::vector<typename Engine::ComplexV> D1_complex_simd(nmax + 1);
  evalDownwardD1<double, Engine>(z_simd, D1_complex_simd);
  
  int target_n = 1;
  Hb_simd = D1_complex_simd[target_n];
  
  auto bn_simd = calc_bn<double, Engine>(target_n, XL_simd, Hb_simd, mL_simd, 
                                         Psi_real_simd[target_n], Zeta_real_simd[target_n], 
                                         Psi_real_simd[target_n-1], Zeta_real_simd[target_n-1]);

  std::vector<double> res_re(lanes), res_im(lanes);
  hn::Store(bn_simd.re, d, res_re.data());
  hn::Store(bn_simd.im, d, res_im.data());

  for (size_t i = 0; i < lanes; ++i) {
    std::complex<double> z_lane(x_vals[i] * mr_vals[i], x_vals[i] * mi_vals[i]);
    std::complex<double> z_real_lane(x_vals[i], 0.0);
    std::complex<double> m_lane(mr_vals[i], mi_vals[i]);
    double XL_lane = x_vals[i];
    
    std::vector<std::complex<double>> D1_real_ref(nmax + 1);
    std::vector<std::complex<double>> D3_real_ref(nmax + 1);
    std::vector<std::complex<double>> PsiZeta_real_ref(nmax + 1);
    std::vector<std::complex<double>> Psi_real_ref(nmax + 1);
    std::vector<std::complex<double>> Zeta_real_ref(nmax + 1);

    evalDownwardD1<double, ScalarEngine>(z_real_lane, D1_real_ref);
    evalUpwardPsi<double, ScalarEngine>(z_real_lane, D1_real_ref, Psi_real_ref);
    evalUpwardD3<double, ScalarEngine>(z_real_lane, D1_real_ref, D3_real_ref, PsiZeta_real_ref);
    for(int k=0; k<=nmax; ++k) Zeta_real_ref[k] = PsiZeta_real_ref[k] / Psi_real_ref[k];

    std::vector<std::complex<double>> D1_complex_ref(nmax + 1);
    evalDownwardD1<double, ScalarEngine>(z_lane, D1_complex_ref);
    std::complex<double> Hb_lane = D1_complex_ref[target_n];

    auto bn_ref = calc_bn<double, ScalarEngine>(target_n, XL_lane, Hb_lane, m_lane,
                                                Psi_real_ref[target_n], Zeta_real_ref[target_n],
                                                Psi_real_ref[target_n-1], Zeta_real_ref[target_n-1]);

    EXPECT_NEAR(res_re[i], bn_ref.real(), 1e-12) << "Lane " << i << " mismatch at bn[" << target_n << "]";
    EXPECT_NEAR(res_im[i], bn_ref.imag(), 1e-12) << "Lane " << i << " mismatch at bn[" << target_n << "]";
  }
}

} // namespace HWY_NAMESPACE
} // namespace nmie

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


