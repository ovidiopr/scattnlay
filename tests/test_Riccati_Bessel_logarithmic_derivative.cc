#include "gtest/gtest.h"
#include "../src/nmie-basic.hpp"
#include "test_spec_functions_data.hpp"
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



void parse_mpmath_data(const double min_abs_tol, const std::tuple< std::complex<double>, int, std::complex<double>, double, double > data,
                       std::complex<double> &z, unsigned int &n, std::complex<double> &func_mp,
                       double &re_abs_tol, double &im_abs_tol){
  z = std::get<0>(data);
  n = std::get<1>(data);
  func_mp = std::get<2>(data);
  re_abs_tol = ( std::get<3>(data) > min_abs_tol && std::real(func_mp) < min_abs_tol)
                    ? std::get<3>(data) : min_abs_tol;
  im_abs_tol = ( std::get<4>(data) > min_abs_tol && std::imag(func_mp) < min_abs_tol)
                    ? std::get<4>(data) : min_abs_tol;
  // if re(func_mp) < 0.5 then round will give 0. To avoid zero tolerance add one.
  re_abs_tol *= std::abs(std::round(std::real(func_mp))) + 1;
  im_abs_tol *= std::abs(std::round(std::imag(func_mp))) + 1;
}
void parse2_mpmath_data(const nmie::FloatType min_abs_tol,
                        const std::tuple< nmie::FloatType, std::complex<nmie::FloatType>, int, std::complex<nmie::FloatType>, nmie::FloatType, nmie::FloatType > data,
                        nmie::FloatType &x, std::complex<nmie::FloatType> &m, unsigned int &n, std::complex<nmie::FloatType> &func_mp,
                        nmie::FloatType &re_abs_tol, nmie::FloatType &im_abs_tol){
  x = std::get<0>(data);
  m = std::get<1>(data);
  n = std::get<2>(data);
  func_mp = std::get<3>(data);
  re_abs_tol = ( std::get<4>(data) > min_abs_tol && std::real(func_mp) < min_abs_tol)
               ? std::get<4>(data) : min_abs_tol;
  im_abs_tol = ( std::get<5>(data) > min_abs_tol && std::imag(func_mp) < min_abs_tol)
               ? std::get<5>(data) : min_abs_tol;
  // if re(func_mp) < 0.5 then round will give 0. To avoid zero tolerance add one.
  re_abs_tol *= std::abs(std::round(std::real(func_mp))) + 1;
  im_abs_tol *= std::abs(std::round(std::imag(func_mp))) + 1;
}

template<class T> inline T pow2(const T value) {return value*value;}

//TEST(an_test, DISABLED_mpmath_generated_input) {
TEST(an_test, mpmath_generated_input) {
  double min_abs_tol = 5e-14, x;
  std::complex<double> m, an_mp;
  unsigned int n;
  double re_abs_tol,  im_abs_tol;
  for (const auto &data : an_test_30digits) {
    parse2_mpmath_data(min_abs_tol, data, x, m, n, an_mp, re_abs_tol, im_abs_tol);
    auto Nstop = nmie::LeRu_near_field_cutoff(m * x)+1;

    nmie::MultiLayerMie<nmie::FloatType> ml_mie;
    ml_mie.SetLayersSize({x});
    ml_mie.SetLayersIndex({m});
    ml_mie.SetMaxTerms(Nstop);
    ml_mie.calcScattCoeffs();
    auto an = ml_mie.GetAn();
//    auto bn = ml_mie.GetBn();

    if (n > an.size()) continue;
    if (n == 0) continue;
    EXPECT_NEAR(std::real(an[n-1]), std::real(an_mp), re_abs_tol)
              << "Db at n=" << n << " Nstop="<< Nstop<<" m="<<m<<" x="<<x;
    EXPECT_NEAR(std::imag(an[n-1]), std::imag(an_mp), im_abs_tol)
              << "Db at n=" << n << " Nstop="<< Nstop<<" m="<<m<<" x="<<x;
  }
}

//TEST(bn_test, DISABLED_mpmath_generated_input) {
TEST(bn_test, mpmath_generated_input) {
  double min_abs_tol = 3e-14, x;
  std::complex<double> m, bn_mp;
  unsigned int n;
  double re_abs_tol,  im_abs_tol;
  for (const auto &data : bn_test_30digits) {
    parse2_mpmath_data(min_abs_tol, data, x, m, n, bn_mp, re_abs_tol, im_abs_tol);
    auto Nstop = nmie::LeRu_near_field_cutoff(m * x)+1;

    nmie::MultiLayerMie<nmie::FloatType> ml_mie;
    ml_mie.SetLayersSize({x});
    ml_mie.SetLayersIndex({m});
    ml_mie.SetMaxTerms(Nstop);
    ml_mie.calcScattCoeffs();
//    auto an = ml_mie.GetAn();
    auto bn = ml_mie.GetBn();

    if (n > bn.size()) continue;
    if (n == 0) continue;
    EXPECT_NEAR(std::real(bn[n-1]), std::real(bn_mp), re_abs_tol)
              << "Db at n=" << n << " Nstop="<< Nstop<<" m="<<m<<" x="<<x;
    EXPECT_NEAR(std::imag(bn[n-1]), std::imag(bn_mp), im_abs_tol)
              << "Db at n=" << n << " Nstop="<< Nstop<<" m="<<m<<" x="<<x;
  }
}

//TEST(zeta_psizeta_test, DISABLED_mpmath_generated_input) {
TEST(zeta_psizeta_test, mpmath_generated_input) {
  double min_abs_tol = 2e-10;
  std::complex<double> z, zeta_mp;
  unsigned int n;
  double re_abs_tol,  im_abs_tol;
  for (const auto &data : zeta_test_16digits) {
    parse_mpmath_data(min_abs_tol, data, z, n, zeta_mp, re_abs_tol, im_abs_tol);
    auto Nstop = nmie::LeRu_near_field_cutoff(z)+10000;
    if (n > Nstop) continue;
    std::vector<std::complex<nmie::FloatType>> D1dr(Nstop+135), D3(Nstop+135),
        PsiZeta(Nstop+135), Psi(Nstop);
    nmie::evalDownwardD1<nmie::FloatType>(z, D1dr);
    nmie::evalUpwardD3<nmie::FloatType>(z, D1dr, D3, PsiZeta);
    nmie::evalUpwardPsi<nmie::FloatType>(z, D1dr, Psi);
    auto a = std::real(PsiZeta[n]);
    auto b = std::imag(PsiZeta[n]);
    auto c = std::real(Psi[n]);
    auto d = std::imag(Psi[n]);
    auto c_one = std::complex<nmie::FloatType>(0, 1);
    auto zeta = (a*c + b*d)/(pow2(c) + pow2(d)) +
        c_one * ((b*c - a*d)/(pow2(c) + pow2(d)));
//    zeta = PsiZeta[n]/Psi[n];
    if (std::isnan(std::real(zeta)) || std::isnan(std::imag(zeta))) continue;
//    std::vector<std::complex<nmie::FloatType>> D1dr(Nstop+35), D3(Nstop+35), zeta(Nstop);
//    nmie::evalDownwardD1<nmie::FloatType>(z, D1dr);
//    nmie::evalUpwardD3(z, D1dr, D3);
//    nmie::evalUpwardZeta(z, D3, zeta);

    EXPECT_NEAR(std::real(zeta), std::real(zeta_mp), re_abs_tol)
              << "zeta at n=" << n << " Nstop="<< Nstop<<" z="<<z;
    EXPECT_NEAR(std::imag(zeta), std::imag(zeta_mp), im_abs_tol)
              << "zeta at n=" << n << " Nstop="<< Nstop<<" z="<<z;
  }
}

// // Old way to evaluate Zeta
// TEST(zeta_test, DISABLED_mpmath_generated_input) {
// //TEST(zeta_test, mpmath_generated_input) {
//  double min_abs_tol = 2e-5;
//  std::complex<double> z, zeta_mp;
//  int n;
//  double re_abs_tol,  im_abs_tol;
//  for (const auto &data : zeta_test_16digits) {
//    parse_mpmath_data(min_abs_tol, data, z, n, zeta_mp, re_abs_tol, im_abs_tol);
//    auto Nstop = nmie::LeRu_near_field_cutoff(z)+10000;
//    if (n > Nstop) continue;
//    std::vector<std::complex<nmie::FloatType>> D1dr(Nstop), D3(Nstop),
//        PsiZeta(Nstop), zeta(Nstop);
//    nmie::evalDownwardD1<nmie::FloatType>(z, D1dr);
//    nmie::evalUpwardD3<nmie::FloatType>(z, D1dr, D3, PsiZeta);
//    nmie::evalUpwardZeta(z, D3, zeta);
//    if (std::isnan(std::real(zeta[n])) || std::isnan(std::imag(zeta[n]))) continue;
//
//    EXPECT_NEAR(std::real(zeta[n]), std::real(zeta_mp), re_abs_tol)
//              << "zeta[n] at n=" << n << " Nstop="<< Nstop<<" z="<<z;
//    EXPECT_NEAR(std::imag(zeta[n]), std::imag(zeta_mp), im_abs_tol)
//              << "zeta at n=" << n << " Nstop="<< Nstop<<" z="<<z;
//  }
//}


//TEST(psizeta_test, DISABLED_mpmath_generated_input) {
TEST(psizeta_test, mpmath_generated_input) {
  double min_abs_tol = 9e-11;
  std::complex<double> z, PsiZeta_mp;
  unsigned int n;
  double re_abs_tol,  im_abs_tol;
  for (const auto &data : psi_mul_zeta_test_16digits) {
    parse_mpmath_data(min_abs_tol, data, z, n, PsiZeta_mp, re_abs_tol, im_abs_tol);
    auto Nstop = nmie::LeRu_near_field_cutoff(z)+10000;
    if (n > Nstop) continue;
    std::vector<std::complex<nmie::FloatType>> D1dr(Nstop), D3(Nstop), PsiZeta(Nstop);
    nmie::evalDownwardD1<nmie::FloatType>(z, D1dr);
    nmie::evalUpwardD3<nmie::FloatType>(z, D1dr, D3, PsiZeta);

    EXPECT_NEAR(std::real(PsiZeta[n]), std::real(PsiZeta_mp), re_abs_tol)
              << "PsiZeta at n=" << n << " Nstop="<< Nstop<<" z="<<z;
    EXPECT_NEAR(std::imag(PsiZeta[n]), std::imag(PsiZeta_mp), im_abs_tol)
              << "PsiZeta at n=" << n << " Nstop="<< Nstop<<" z="<<z;
//    std::vector<nmie::FloatType> PsiUp(Nstop);
//    nmie::evalPsi(std::real(z), PsiUp);
//    EXPECT_NEAR(((PsiUp[n])), std::real(PsiZeta_mp), re_abs_tol)
//              << "PsiZeta(up) at n=" << n << " z="<<z;
  }
}


TEST(psi_test, mpmath_generated_input) {
  double min_abs_tol = 1e-12;
  std::complex<double> z, Psi_mp;
  unsigned int n;
  double re_abs_tol,  im_abs_tol;
  for (const auto &data : psi_test_16digits) {
    parse_mpmath_data(min_abs_tol, data, z, n, Psi_mp, re_abs_tol, im_abs_tol);
    auto Nstop = nmie::LeRu_near_field_cutoff(z)+10000;
    if (n > Nstop) continue;
    std::vector<std::complex<nmie::FloatType>> D1dr(Nstop+35), Psi(Nstop);
    nmie::evalDownwardD1<nmie::FloatType>(z, D1dr);
    nmie::evalUpwardPsi<nmie::FloatType>(z, D1dr, Psi);

    EXPECT_NEAR(std::real(Psi[n]), std::real(Psi_mp), re_abs_tol)
              << "Psi at n=" << n << " Nstop="<< Nstop<<" z="<<z;
    EXPECT_NEAR(std::imag(Psi[n]), std::imag(Psi_mp), im_abs_tol)
              << "Psi at n=" << n << " Nstop="<< Nstop<<" z="<<z;
  }
}


//TEST(D3test, DISABLED_mpmath_generated_input) {
TEST(D3test, mpmath_generated_input) {
  double min_abs_tol = 2e-11;
  std::complex<double> z, D3_mp;
  unsigned int n;
  double re_abs_tol,  im_abs_tol;
  for (const auto &data : D3_test_16digits) {
    parse_mpmath_data(min_abs_tol, data, z, n, D3_mp, re_abs_tol, im_abs_tol);
    auto Nstop = nmie::LeRu_near_field_cutoff(z)+35;
    std::vector<std::complex<nmie::FloatType>> D1dr(Nstop), D3(Nstop), PsiZeta(Nstop);
    nmie::evalDownwardD1<nmie::FloatType>(z, D1dr);
    nmie::evalUpwardD3<nmie::FloatType>(z, D1dr, D3, PsiZeta);

    EXPECT_NEAR(std::real(D3[n]), std::real(D3_mp), re_abs_tol)
              << "D3 at n=" << n << " Nstop="<< Nstop<<" z="<<z;
    EXPECT_NEAR(std::imag(D3[n]), std::imag(D3_mp), im_abs_tol)
              << "D3 at n=" << n << " Nstop="<< Nstop<<" z="<<z;

  }
}


//TEST(D1test, DISABLED_mpmath_generated_input) {
  TEST(D1test, mpmath_generated_input) {
  double min_abs_tol = 7e-11, x;
  std::complex<double> m, z, D1_mp;
  unsigned int n;
  double re_abs_tol,  im_abs_tol;
  for (const auto &data : D1_test_30digits) {
    parse2_mpmath_data(min_abs_tol, data, x, m, n, D1_mp, re_abs_tol, im_abs_tol);
    if (n == 0 && nmie::cabs(D1_mp) > 1e14) continue;
    z = m*x;
//    auto Nstop = nmie::LeRu_near_field_cutoff(z)+1;
//    auto Nstop = n;
    int valid_digits = 16;
    int nstar = nmie::getNStar<nmie::FloatType>(n, z, valid_digits);
    std::vector<std::complex<nmie::FloatType>> Db(nstar),Dold(nstar), r;
    r.resize(nstar);
    Db.resize(nstar);
    nmie::evalBackwardR(z,r);
    nmie::convertRtoD1(z, r, Db);
    if (n > Db.size()) continue;
    EXPECT_NEAR(std::real(Db[n]), std::real(D1_mp), re_abs_tol)
              << "Db at n=" << n <<" nstar="<<nstar<< " z="<<z;
    EXPECT_NEAR(std::imag(Db[n]), std::imag(D1_mp), im_abs_tol)
              << "Db at n=" << n <<" nstar="<<nstar<< " z="<<z;
    nmie::evalDownwardD1<nmie::FloatType>(z, Dold);
    if (n > Dold.size()) continue;
    EXPECT_NEAR(std::real(Dold[n]), std::real(D1_mp), re_abs_tol)
              << "Dold at n=" << n << " z="<<z;
    EXPECT_NEAR(std::imag(Dold[n]), std::imag(D1_mp), im_abs_tol)
              << "Dold at n=" << n << " z="<<z;

  }
}


//TEST(D1test, DISABLED_WYang_data){
TEST(D1test, WYang_data){
  double abs_tol = 4e-10;
  int test_loss_digits = std::round(15 - std::log10(1/abs_tol));
  int Nstop = 131;
  std::vector<std::complex<nmie::FloatType>> Df(Nstop), Db(Nstop),Dold(Nstop), r;
  std::complex<nmie::FloatType> z(1.05,1);
  z = z*80.0;
// eval D1 directly from backward recurrence
  nmie::evalDownwardD1<nmie::FloatType>(z, Dold);

//  eval forward recurrence
  r.resize(Nstop+1);
  nmie::evalForwardR(z, r);
  nmie::convertRtoD1(z, r, Df);

for (unsigned int i = 0; i < Dtest_n.size(); i++) {
  unsigned int n = Dtest_n[i];
  int forward_loss_digits = nmie::evalKapteynNumberOfLostSignificantDigits<nmie::FloatType>(n, z);
  forward_loss_digits += 3; // Kapteyn is too optimistic
  if (test_loss_digits > forward_loss_digits ) {
    EXPECT_NEAR(std::real(Df[n]), std::real(Dtest_D1[i]),
                abs_tol) << "f at n=" << n << " lost digits = " << forward_loss_digits;
    EXPECT_NEAR(std::imag(Df[n]), std::imag(Dtest_D1[i]),
                abs_tol) << "f at n=" << n << " lost digits = " << forward_loss_digits;
  }
// eval backward recurrence
  int valid_digits = 6;
  int nstar = nmie::getNStar<nmie::FloatType>(n, z, valid_digits);
  r.resize(nstar);
  Db.resize(nstar);
  nmie::evalBackwardR(z,r);
  nmie::convertRtoD1(z, r, Db);

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
  double l = nmie::evalKapteynNumberOfLostSignificantDigits<double>(80, std::complex<double>(100,100));
  EXPECT_EQ(l, 7)<<"Should be equal";
  std::complex<double> z(10000,0);
  l = nmie::evalKapteynNumberOfLostSignificantDigits<double>(5070, z);
  EXPECT_EQ(l, 0)<<"Should be equal";
  // find NStar such that l_nstar(z) - l_nmax(z) >= valid_digits
  int NStar = nmie::getNStar<double>(5070, z,6);
  EXPECT_GE(NStar, 10130);
//  const double pi=3.14159265358979323846;
//  z = std::complex<double>(100,100);
//  l = nmie::evalKapteynNumberOfLostSignificantDigits(1, z);
//  EXPECT_EQ(l, 0)<<"Should be equal";

}


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
