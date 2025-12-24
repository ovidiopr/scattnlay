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

} // namespace HWY_NAMESPACE
} // namespace nmie

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


