#include "gtest/gtest.h"
#include "../src/nmie-impl.hpp"
// From W. Yang APPLIED OPTICS Vol. 42, No. 9,  20 March 2003
// Dtest refractive index is m={1.05,1}, the size parameter is x = 80
std::vector<int> Dtest_n({0,1,30,50,60,70,75,80,85,90,99,116,130});
std::vector< std::complex<double>>
    Dtest_D1({{0.11449e-15 ,-0.10000e+01 },{0.74646e-04 ,-0.10000e+01 },
              {0.34764e-01 ,-0.99870},{0.95292e-01 ,-0.99935},
              {0.13645,-0.10019e+01 },{0.18439,-0.10070e+01 },
              {0.21070,-0.10107e+01 },{0.23845,-0.10154e+01 },
              {0.26752,-0.10210e+01 },{0.29777,-0.10278e+01 },
              {0.35481,-0.10426e+01 },{0.46923,-0.10806e+01 },
              {0.17656,-0.13895e+01 }});
std::vector< std::complex<double>>
    Dtest_D2({{0.64966e-69 ,-0.10000e+01 },{0.74646e-04 ,-0.10000e+01 },
              {0.34764e-01 ,-0.99870},{0.95292e-01 ,-0.99935},
              {0.13645,-0.10019e+01 },{0.17769,-0.10099e+01 },
              {0.41264e-01 ,-0.21076e+01 },{-0.20190,0.10435e+01 },
              {-0.26343,0.10223e+01 },{-0.29339,0.10291e+01 },
              {-0.34969,0.10437e+01 },{-0.46296,0.10809e+01 },
              {-0.56047,0.11206e+01 }});
std::vector< std::complex<double>>
    Dtest_D3({{0.00000,0.10000e+01 },{-0.73809e-04 ,0.10000e+01 },
              {-0.34344e-01 ,0.99912},{-0.94022e-01 ,0.10004e+01 },
              {-0.13455,0.10032e+01 },{-0.18172,0.10084e+01 },
              {-0.20762,0.10122e+01 },{-0.23494,0.10169e+01 },
              {-0.26357,0.10225e+01 },{-0.29339,0.10291e+01 },
              {-0.34969,0.10437e+01 },{-0.46296,0.10809e+01 },
              {-0.56047,0.11206e+01 }});


//TEST(Dtest, DISABLED_WYang_data){
TEST(Dtest, WYang_data){
  double rel_tol = 1e-6;
  int Nstop = 131;
  std::vector<std::complex<nmie::FloatType>> Df(Nstop), Db(Nstop), r;
  std::complex<nmie::FloatType> z(1.05,1);
  z = z*80.0;
  nmie::evalForwardD1(z, Df);
  r.resize(Nstop+1);
  nmie::evalForwardR(z, r);
  nmie::convertRtoD1(z, r, Df);
  int valid_digits = 6;
  int nstar = nmie::getNStar(Nstop, z, valid_digits);
  r.resize(nstar);
  nmie::evalBackwardR(z,r);
  nmie::convertRtoD1(z, r, Db);

for (int i = 0; i < Dtest_n.size(); i++) {
  int n = Dtest_n[i];
//  EXPECT_FLOAT_EQ(std::real(Df[n]), std::real(Dtest_D1[i])) << "f at n=" << n;
//  EXPECT_FLOAT_EQ(std::imag(Df[n]), std::imag(Dtest_D1[i])) << "f at n=" << n;
  EXPECT_NEAR(std::real(Db[n]), std::real(Dtest_D1[i]),
              rel_tol*std::real(Db[n])) << "b at n=" << n;
  EXPECT_NEAR(std::imag(Db[n]), std::imag(Dtest_D1[i]),
              rel_tol*std::imag(Db[n])) << "b at n=" << n;
//  EXPECT_FLOAT_EQ(std::real(Df[n]), std::real(Db[n])) << "f-b at n=" << n;
//  EXPECT_FLOAT_EQ(std::imag(Df[n]), std::imag(Db[n])) << "f-b at n=" << n;
  }
}

TEST(ratio_funcion_forward_vs_backward_recurrence, compare){
  int Nstop = 131;
  std::vector<std::complex<nmie::FloatType>> rf(Nstop), rb;
//  std::complex<nmie::FloatType> z(1.05,1);
//  z = z*80.0;

  std::complex<nmie::FloatType> z(1.3,-2.1);
// rb[49] in Wolfram Alpha
// n=49, z=1.3-2.1i,  SphericalBesselJ[n-1,z]/SphericalBesselJ[n,z]
  nmie::evalForwardR(z, rf);
  int valid_digits = 20;
  int nstar = nmie::getNStar(Nstop, z, valid_digits);
  rb.resize(nstar);
  nmie::evalBackwardR(z,rb);

  for (int i = 0; i < Dtest_n.size(); i++) {
    int n = Dtest_n[i];
    EXPECT_FLOAT_EQ(std::real(rf[n]), std::real(rb[i]))
              << "at n=" << n
              << " expected forward loss="<<nmie::evalKapteynNumberOfLostSignificantDigits(n,z);
    EXPECT_FLOAT_EQ(std::imag(rf[n]), std::imag(rb[i]))
              << "at n=" << n
              << " expected forward loss="<<nmie::evalKapteynNumberOfLostSignificantDigits(n,z);
  }
}


TEST(KaptyenTest, HandlesInput) {
  // H.Du APPLIED OPTICS, Vol. 43, No. 9, 20 March 2004
  double l = nmie::evalKapteynNumberOfLostSignificantDigits(80, std::complex<double>(100,100));
  EXPECT_EQ(l, 7)<<"Should be equal";
  std::complex<double> z(10000,0);
  l = nmie::evalKapteynNumberOfLostSignificantDigits(5070, z);
  EXPECT_EQ(l, 0)<<"Should be equal";
  int NStar = nmie::getNStar(5070, z,6);
  EXPECT_GE(NStar, 10130);
}


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
