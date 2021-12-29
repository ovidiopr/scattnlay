#include "gtest/gtest.h"
#include "../src/nmie-basic.hpp"
#include "../src/nmie-nearfield.hpp"
#include "test_cases.hpp"

//TEST(RunFieldCalculationCartesian, DISABLED_HandlesInput) {
TEST(RunFieldCalculationCartesian, HandlesInput) {
  nmie::MultiLayerMie<nmie::FloatType> nmie;
//  EXPECT_THROW(nmie.RunFieldCalculationPolar(0), std::invalid_argument);
//  EXPECT_THROW(nmie.RunFieldCalculationPolar(1,1,10,5), std::invalid_argument);
  double r = 2*nmie.PI_*100/619;
//  double r = 1500;
  nmie.SetLayersSize({r/2, r});
  nmie.SetLayersIndex({ {4.0,0}, {4.0,0}});
  nmie.RunMieCalculation();
  int nmax = 21;
  // TODO add check of E and H symmetry for X and Y axis inversion
  nmie.RunFieldCalculationCartesian(2, 2, nmie::Planes::kEk,
                                    0, 0, 0, true, nmax);
  nmie.RunFieldCalculationCartesian(2, 2, nmie::Planes::kHk,
                                    0, 0, 0, true, nmax);
  nmie.RunFieldCalculationCartesian(2, 2, nmie::Planes::kEH,
                                    0, 0, 0, true, nmax);

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
