#include "gtest/gtest.h"
#include "../src/nmie-impl.hpp"
#include "test_spec_functions_data.h"
// From W. Yang APPLIED OPTICS Vol. 42, No. 9,  20 March 2003
// Dtest refractive index is m={1.05,1}, the size parameter is x = 80
std::vector<int> Dtest_n({0,1,30,50,60,70,75,80,85,90,99,116,130});
std::vector< std::complex<double>>
    Dtest_D1({
      //Orig
//                 {0.11449e-15 ,-0.10000e+01 },{0.74646e-04 ,-0.10000e+01 },
//                 {0.34764e-01 ,-0.99870},{0.95292e-01 ,-0.99935},
//                 {0.13645,-0.10019e+01 },{0.18439,-0.10070e+01 },
//                 {0.21070,-0.10107e+01 },{0.23845,-0.10154e+01 },
//                 {0.26752,-0.10210e+01 },{0.29777,-0.10278e+01 },
//                 {0.35481,-0.10426e+01 },{0.46923,-0.10806e+01 },
//                 {0.17656,-0.13895e+01 }
      // mod (from Python mpmath)
                 {0.0,-1.0}, {7.464603828e-5,-0.9999958865},
                 {0.03476380918,-0.9986960672},{0.09529213152,-0.999347654},
                 {0.1364513887,-1.001895883},{0.184388335,-1.006979164},
                 {0.2107044267,-1.01072099},{0.2384524295,-1.015382914},
                 {0.2675164524,-1.021040337},{0.2977711192,-1.027753418},
                 {0.3548096904,-1.042622957},{0.4692294405,-1.080629479},
                 {0.5673827836,-1.121108944},
             });

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


int LeRu_cutoff(std::complex<double> z) {
  auto x = std::abs(z);
  return std::round(x + 11 * std::pow(x, (1.0 / 3.0)) + 1);
}


TEST(D1test, mpmath_generated_input) {
  double min_abs_tol = 1e-13;
  for (const auto &data : D1_test_10digits) {
    auto z = std::get<0>(data);
    auto n = std::get<1>(data);
    auto D1_mpmath = std::get<2>(data);
    auto re_abs_tol = (std::get<3>(data)>min_abs_tol) ? std::get<3>(data) : min_abs_tol;
    auto im_abs_tol = (std::get<4>(data)>min_abs_tol) ? std::get<4>(data) : min_abs_tol;
    // if re(D1) < 0.5 then round will give 0. To avoid zero tolerance add one.
    // To
    re_abs_tol *= std::abs(std::round(std::real(D1_mpmath))) + 11;
    im_abs_tol *= std::abs(std::round(std::imag(D1_mpmath))) + 11;
    auto Nstop = LeRu_cutoff(z)+1;
    std::vector<std::complex<nmie::FloatType>> Df(Nstop), Db(Nstop),Dold(Nstop), r;
    int valid_digits = 6;
    int nstar = nmie::getNStar(Nstop, z, valid_digits);
    r.resize(nstar);
    nmie::evalBackwardR(z,r);
    nmie::convertRtoD1(z, r, Db);
    if (n > Db.size()) continue;
    EXPECT_NEAR(std::real(Db[n]), std::real(D1_mpmath), re_abs_tol)
              << "b at n=" << n << " Nstop="<< Nstop<<" nstar="<<nstar<< " z="<<z;
    EXPECT_NEAR(std::imag(Db[n]), std::imag(D1_mpmath), im_abs_tol)
              << "b at n=" << n << " Nstop="<< Nstop<<" nstar="<<nstar<< " z="<<z;
  }
}


//TEST(D1test, DISABLED_WYang_data){
TEST(D1test, WYang_data){
  double abs_tol = 1e-9;
  int test_loss_digits = std::round(15 - std::log10(1/abs_tol));
  int Nstop = 131;
  std::vector<std::complex<nmie::FloatType>> Df(Nstop), Db(Nstop),Dold(Nstop), r;
  std::complex<nmie::FloatType> z(1.05,1);
  z = z*80.0;
// eval D1 directly from backward recurrence
  nmie::evalDownwardD1(z, Dold);

//  eval forward recurrence
  r.resize(Nstop+1);
  nmie::evalForwardR(z, r);
  nmie::convertRtoD1(z, r, Df);
// eval backward recurrence
  int valid_digits = 6;
  int nstar = nmie::getNStar(Nstop, z, valid_digits);
  r.resize(nstar);
  nmie::evalBackwardR(z,r);
  nmie::convertRtoD1(z, r, Db);

for (int i = 0; i < Dtest_n.size(); i++) {
  int n = Dtest_n[i];
  int forward_loss_digits = nmie::evalKapteynNumberOfLostSignificantDigits(n, z);
  forward_loss_digits += 3; // Kapteyn is too optimistic
  if (test_loss_digits > forward_loss_digits ) {
    EXPECT_NEAR(std::real(Df[n]), std::real(Dtest_D1[i]),
                abs_tol) << "f at n=" << n << " lost digits = " << forward_loss_digits;
    EXPECT_NEAR(std::imag(Df[n]), std::imag(Dtest_D1[i]),
                abs_tol) << "f at n=" << n << " lost digits = " << forward_loss_digits;
  }
  EXPECT_NEAR(std::real(Db[n]), std::real(Dtest_D1[i]),
              abs_tol) << "b at n=" << n;
  EXPECT_NEAR(std::imag(Db[n]), std::imag(Dtest_D1[i]),
              abs_tol) << "b at n=" << n;
  if (n < Dold.size()-15) {
    EXPECT_NEAR(std::real(Dold[n]), std::real(Dtest_D1[i]),
                abs_tol) << "old at n=" << n;
    EXPECT_NEAR(std::imag(Dold[n]), std::imag(Dtest_D1[i]),
                abs_tol) << "old at n=" << n;
  }

}
}



TEST(KaptyenTest, HandlesInput) {
  // H.Du APPLIED OPTICS, Vol. 43, No. 9, 20 March 2004
  double l = nmie::evalKapteynNumberOfLostSignificantDigits(80, std::complex<double>(100,100));
  EXPECT_EQ(l, 7)<<"Should be equal";
  std::complex<double> z(10000,0);
  l = nmie::evalKapteynNumberOfLostSignificantDigits(5070, z);
  EXPECT_EQ(l, 0)<<"Should be equal";
  // find NStar such that l_nstar(z) - l_nmax(z) >= valid_digits
  int NStar = nmie::getNStar(5070, z,6);
  EXPECT_GE(NStar, 10130);
}


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
