#include "gtest/gtest.h"
#include "../src/nmie-basic.hpp"
#include "../src/nmie-nearfield.hpp"
#include "test_cases.hpp"

//TEST(RunFieldCalculationCartesian, DISABLED_HandlesInput) {
TEST(RunFieldCalculationCartesian, HandlesInput) {
  nmie::MultiLayerMie<nmie::FloatType> nmie;
//  EXPECT_THROW(nmie.RunFieldCalculationPolar(0), std::invalid_argument);
//  EXPECT_THROW(nmie.RunFieldCalculationPolar(1,1,10,5), std::invalid_argument);
  nmie::FloatType total_r = 2*nmie.PI_*1000/532;
//  double r = 1500;
  nmie.SetLayersSize({total_r/2, total_r});
  nmie.SetLayersIndex({ {1.330,0}, {1.33,0}});
  nmie.RunMieCalculation();
  double relative_max_distance = 1e-10;
//  nmie.SetModeNmaxAndType(3,-1);
//  int nmax = 21;
  nmie.RunFieldCalculationCartesian(2, 5, relative_max_distance, nmie::Planes::kEk,
                                    1.0, 0, 0, false,3);
  auto Eabs = nmie.GetFieldEabs();
  auto E = nmie.GetFieldE();
  std::cout<<std::endl;
  {
    // Eabs points are located near the sphere outer border
    //
    //    0   1   2   3   4
    //    ----- border ----
    //    5   6   7   8   9
    // distance between points (0) and (4) is relative_max_distance*total_r, initial
    // value used for the test was 1e-10*total_r, so we expect good linear dependence
    // for points from 0 to 4 and 5 to 9. In the asserts we check, that the slope doesn't
    // change too fast inside the curve. While writing this, the test was failing.
    // The value of z-coordinates of 2 and 7 points = 0
    using nmie::nmm::abs;
    EXPECT_TRUE(
        ( abs(Eabs[0] - Eabs[1]) + abs(Eabs[3] - Eabs[4]) ) >= abs(Eabs[1] - Eabs[2])
        );
    EXPECT_TRUE(
        ( abs(Eabs[5] - Eabs[6]) + abs(Eabs[8] - Eabs[9]) ) >= abs(Eabs[6] - Eabs[7])
    );
  }


//  nmie.RunFieldCalculationCartesian(2, 2, 2, nmie::Planes::kHk,
//                                    0, 0, 0, true);
//  nmie.RunFieldCalculationCartesian(2, 2, 2, nmie::Planes::kEH,
//                                    0, 0, 0, true);
  // TODO add check of E and H symmetry for X and Y axis inversion

//  EXPECT_EQ(1, nmie.GetMaxTerms());
//  EXPECT_FALSE(nmie.GetFieldConvergence());
//  auto Eabs = nmie.GetFieldEabs();
//  EXPECT_TRUE(std::isnan(static_cast<double>(Eabs[0])));

}

//TEST(RunFieldCalculationPolar, DISABLED_HandlesInput) {
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
//#ifndef MULTI_PRECISION
//TEST(BulkSphere, DISABLED_HandlesInput) {
TEST(BulkSphere, HandlesInput) {
  nmie::MultiLayerMie<nmie::FloatType> nmie;
  for (const auto &data : parameters_bulk_sphere) {
    auto x = std::get<0>(data);
    auto m = std::get<1>(data);
    nmie.SetLayersSize({x});
    nmie.SetLayersIndex({m});
    nmie.SetMaxTerms(-1);
//    nmie.RunMieCalculation();
//    std::cout<<" test case: "<<std::get<2>(data)<<" Qsca="<<nmie.GetQsca()<<std::endl;
    nmie.RunFieldCalculationPolar(4,3,x,x*3, 0, static_cast<double>(nmie.PI_), 0, static_cast<double>(nmie.PI_),true, -1);
    auto Eabs = nmie.GetFieldEabs();
    for (auto &E:Eabs) E=nmie::pow2(E);
//    print(Eabs)
    EXPECT_TRUE(nmie.GetFieldConvergence())<<"Outside test for x="<<x<<" m="<<m<<" test case: "<<std::get<2>(data)<<std::endl;
    nmie.RunFieldCalculationPolar(4,10,x*0.01,x, 0, static_cast<double>(nmie.PI_), 0, static_cast<double>(nmie.PI_),true, -1);
    EXPECT_TRUE(nmie.GetFieldConvergence())<<"Inside test for x="<<x<<" m="<<m<<" test case: "<<std::get<2>(data)<<std::endl;
  }
}
//#endif

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
