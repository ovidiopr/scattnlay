#include <iostream>
#include <vector>
#include "hwy/highway.h"
#include "../src/nmie-precision.hpp"

// Simplified smoke test using static dispatch only
// to verify Highway build chain and basic SIMD operations.

namespace nmie {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

void SmokeTestSIMD() {
  const hn::ScalableTag<double> d;
  auto v1 = hn::Set(d, 1.5);
  auto v2 = hn::Set(d, 2.5);
  auto res = hn::Add(v1, v2);
  
  std::vector<double> out(hn::Lanes(d));
  hn::Store(res, d, out.data());

  if (out[0] == 4.0) {
    std::cout << "SIMD Smoke Test Passed: 1.5 + 2.5 = " << out[0] 
              << " on target " << hwy::TargetName(HWY_TARGET) << std::endl;
  } else {
    std::cerr << "SIMD Smoke Test FAILED!" << std::endl;
    exit(1);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace nmie

int main() {
  // Call the static target implementation directly
  nmie::HWY_NAMESPACE::SmokeTestSIMD();
  return 0;
}
