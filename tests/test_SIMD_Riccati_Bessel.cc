#include "gtest/gtest.h"
#include "hwy/highway.h"
#include "../src/special-functions-impl.hpp"
#include "../src/nmie-precision.hpp"

// #define HWY_TARGET_INCLUDE "tests/test_SIMD_Riccati_Bessel.cc"
// #include "hwy/foreach_target.h"

namespace nmie {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

TEST(SIMDRiccatiBessel, ComplexCotSmoke) {
  const hn::ScalableTag<double> d;
  using Engine = HighwayEngine<double>;
  
  // Prepare input: z = 1.0 + 0.5i
  auto z_re = hn::Set(d, 1.0);
  auto z_im = hn::Set(d, 0.5);
  typename Engine::ComplexV z_batch{z_re, z_im};

  // Currently complex_cot in special-functions-impl.hpp expects std::complex.
  // We will need to update the signature in the next step to accept Engine::ComplexV.
  // This test will fail compilation until Step 2.3.
  auto result = complex_cot<double, Engine>(z_batch);
  (void)result;
}

} // namespace HWY_NAMESPACE
} // namespace nmie

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


