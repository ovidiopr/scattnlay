#include "gtest/gtest.h"
#include "hwy/highway.h"
#include "../src/special-functions-impl.hpp"
#include "../src/nmie-precision.hpp"

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

} // namespace HWY_NAMESPACE
} // namespace nmie

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


