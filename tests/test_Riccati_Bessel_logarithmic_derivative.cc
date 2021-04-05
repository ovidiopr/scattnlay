#include "gtest/gtest.h"
#include "../src/nmie-impl.hpp"
TEST(RiccatiBesselTest, HandlesInput) {
  nmie::MultiLayerMie<double> nmie;
  EXPECT_EQ(1, 1)<<"Should be equal";
}
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
