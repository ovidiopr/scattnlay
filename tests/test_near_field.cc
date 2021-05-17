#include "gtest/gtest.h"
#include "../src/nmie-basic.hpp"
#include "../src/nmie-nearfield.hpp"
#include "test_cases.hpp"

TEST(RunFieldCalculationPolar, HandlesInput) {
  nmie::MultiLayerMie<nmie::FloatType> nmie;
  EXPECT_THROW(nmie.RunFieldCalculationPolar(0), std::invalid_argument);
  EXPECT_THROW(nmie.RunFieldCalculationPolar(1,1,10,5), std::invalid_argument);
  double r = 60;
//  double r = 1500;
  nmie.SetLayersSize({r/2, r});
  nmie.SetLayersIndex({ {1.33,0}, {1.33,0}});
  nmie.RunMieCalculation();
  nmie.RunFieldCalculationPolar(1, 1,
                                0.5145*r,
                                r*0.5148,
                                0, 3.14, 0, 0, true, 1);
  EXPECT_EQ(1, nmie.GetMaxTerms());
  EXPECT_FALSE(nmie.GetFieldConvergence());
  auto Eabs = nmie.GetFieldEabs();
  EXPECT_TRUE(std::isnan(static_cast<double>(Eabs[0])));

}
#ifndef MULTI_PRECISION
TEST(BulkSphere, HandlesInput) {
  nmie::MultiLayerMie<nmie::FloatType> nmie;
  for (const auto &data : parameters_bulk_sphere) {
    auto x = std::get<0>(data);
    auto m = std::get<1>(data);
    nmie.SetLayersSize({x});
    nmie.SetLayersIndex({m});
    nmie.RunFieldCalculationPolar(4,3,x,x*3);
    EXPECT_TRUE(nmie.GetFieldConvergence())<<"for x="<<x<<" m="<<m<<" test case: "<<std::get<2>(data)<<std::endl;
    nmie.RunFieldCalculationPolar(4,10,x*0.01,x);
    EXPECT_TRUE(nmie.GetFieldConvergence())<<"for x="<<x<<" m="<<m<<" test case: "<<std::get<2>(data)<<std::endl;
  }
}
#endif

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
