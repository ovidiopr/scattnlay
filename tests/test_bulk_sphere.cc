#include <complex>
#include "../src/mesomie.hpp"
#include "../src/nmie-basic.hpp"
#include "gtest/gtest.h"

// TODO fails for MP with 100 digits. And 16 digits, which should be equal to
// double precision.
#ifndef MULTI_PRECISION
// TEST(BulkSphere, DISABLED_ArgPi) {
//******************************************************************************
TEST(BulkSphere, ArgPi) {
  std::vector<double> WLs{50, 80, 100, 200, 400};  // nm
  double host_index = 2.;
  double core_radius = 100.;  // nm
  double delta = 1e-5;
  nmie::MultiLayerMie<nmie::FloatType> nmie;
  nmie.SetLayersIndex({std::complex<double>(4, 0)});
  for (auto WL : WLs) {
    nmie.SetLayersSize(
        {2 * nmie::PI_ * host_index * core_radius / (WL + delta)});
    nmie.RunMieCalculation();
    double Qabs_p = std::abs(static_cast<double>(nmie.GetQabs()));

    nmie.SetLayersSize(
        {2 * nmie::PI_ * host_index * core_radius / (WL - delta)});
    nmie.RunMieCalculation();
    double Qabs_m = std::abs(static_cast<double>(nmie.GetQabs()));

    nmie.SetLayersSize({2 * nmie::PI_ * host_index * core_radius / (WL)});
    nmie.RunMieCalculation();
    double Qabs = std::abs(static_cast<double>(nmie.GetQabs()));
    EXPECT_GT(Qabs_p + Qabs_m, Qabs);
  }
}
#endif

//******************************************************************************
// A list of tests for a bulk sphere from
// Hong Du, "Mie-scattering calculation," Appl. Opt. 43, 1951-1956 (2004)
// table 1: sphere size and refractive index
// followed by resulting extinction and scattering efficiencies
std::vector<std::tuple<double, std::complex<double>, double, double, char> >
    parameters_and_results{
        // x, {Re(m), Im(m)}, Qext, Qsca, test_name
        {0.099, {0.75, 0}, 7.417859e-06, 7.417859e-06, 'a'},
        {0.101, {0.75, 0}, 8.033538e-06, 8.033538e-06, 'b'},
        {10, {0.75, 0}, 2.232265, 2.232265, 'c'},
        {1000, {0.75, 0}, 1.997908, 1.997908, 'd'},
        {100, {1.33, 1e-5}, 2.101321, 2.096594, 'e'},
        {10000, {1.33, 1e-5}, 2.004089, 1.723857, 'f'},
        {0.055, {1.5, 1}, 0.10149104, 1.131687e-05, 'g'},
        {0.056, {1.5, 1}, 0.1033467, 1.216311e-05, 'h'},
        {100, {1.5, 1}, 2.097502, 1.283697, 'i'},
        {10000, {1.5, 1}, 2.004368, 1.236574, 'j'},
        {1, {10, 10}, 2.532993, 2.049405, 'k'},
        {100, {10, 10}, 2.071124, 1.836785, 'l'},
        {10000, {10, 10}, 2.005914, 1.795393, 'm'},
        //  {1000000, {10, 10}, 2.00022, 1.79218, 'x'}
        // Qsca 1.792181 - from Scattnlay 2025-12-28
    };
//******************************************************************************
// Helper functions for significant figure rounding
double RoundToSigFigs(double val, int n) {
  if (val == 0.0)
    return 0.0;
  double multiplier = std::pow(10, n - std::ceil(std::log10(std::abs(val))));
  return std::round(val * multiplier) / multiplier;
}

std::complex<double> RoundComplex(std::complex<double> val, int n) {
  return {RoundToSigFigs(val.real(), n), RoundToSigFigs(val.imag(), n)};
}

struct DuAmplitudeCase {
  double x;
  std::complex<double> m;
  std::complex<double> s1_0;
  std::complex<double> s1_pi;
  char id;
};

#ifndef MULTI_PRECISION
TEST(BulkSphere, DuTable2Amplitudes) {
  std::vector<DuAmplitudeCase> cases = {
      {0.099,
       {0.75, 0.0},
       {1.81756e-8, 1.65423e-4},
       {1.81756e-8, 1.64810e-4},
       'a'},
      {0.101,
       {0.75, 0.0},
       {2.04875e-8, 1.75642e-4},
       {2.04875e-8, 1.74965e-4},
       'b'},
      {10.0, {0.75, 0.0}, {55.8066, 9.75810}, {-1.07857, 0.0360881}, 'c'},
      {1000.0, {0.75, 0.0}, {499477.0, 13365.0}, {17.0578, -484.251}, 'd'},
      {100.0, {1.33, 1e-5}, {5253.3, 124.319}, {-56.5921, -46.5097}, 'e'},
      // {10000.0, {1.33, 1e-5}, {5.01022e7, 153582.0}, {-182.119, 951.912},
      // 'f'}, //TODO explore minor discrepancies
      {0.055,
       {1.5, 1.0},
       {7.67526e-5, -8.34388e-5},
       {7.66140e-5, -8.33814e-5},
       'g'},
      {0.056,
       {1.5, 1.0},
       {8.10238e-5, -8.80725e-5},
       {8.08721e-5, -8.80098e-5},
       'h'},
      {100.0, {1.5, 1.0}, {5243.75, 293.417}, {-20.2936, -4.38444}, 'i'},
      {10000.0, {1.5, 1.0}, {5.01092e7, 175340.0}, {-218.472, 2064.61}, 'j'},
      // {1.0, {10.0, 10.0}, {0.633248, -0.417931}, {0.448546,
      // -0.791236}, 'k'}, //TODO explore minor discrepancies
      {100.0, {10.0, 10.0}, {5177.81, 26.3381}, {-41.4538, 18.2181}, 'l'},
      {10000.0, {10.0, 10.0}, {5.01479e7, 120600.0}, {2252.48, 3924.47}, 'm'}};

  nmie::MultiLayerMie<nmie::FloatType> nmie;
  std::vector<nmie::FloatType> theta = {0.0, nmie::PI_};  // 0 and 180 degrees

  for (const auto& c : cases) {
    nmie.SetLayersSize({static_cast<nmie::FloatType>(c.x)});
    nmie.SetLayersIndex(
        {std::complex<nmie::FloatType>(c.m.real(), c.m.imag())});
    nmie.SetAngles(theta);
    nmie.RunMieCalculation();

    auto S1 = nmie.GetS1();
    auto S2 = nmie.GetS2();

    // Round results to 6 significant figures to match table format
    auto s1_0_val = std::complex<double>(static_cast<double>(S1[0].real()),
                                         static_cast<double>(S1[0].imag()));
    auto s1_pi_val = std::complex<double>(static_cast<double>(S1[1].real()),
                                          static_cast<double>(S1[1].imag()));

    auto calc_s1_0 = RoundComplex(s1_0_val, 6);
    auto calc_s1_pi = RoundComplex(s1_pi_val, 6);

    EXPECT_NEAR(calc_s1_0.real(), c.s1_0.real(), 1e-15)
        << "Fail Case " << c.id << " S1(0) real";
    EXPECT_NEAR(calc_s1_0.imag(), c.s1_0.imag(), 1e-15)
        << "Fail Case " << c.id << " S1(0) imag";

    EXPECT_NEAR(calc_s1_pi.real(), c.s1_pi.real(), 1e-15)
        << "Fail Case " << c.id << " S1(pi) real";
    EXPECT_NEAR(calc_s1_pi.imag(), c.s1_pi.imag(), 1e-15)
        << "Fail Case " << c.id << " S1(pi) imag";

    // S1(0) = S2(0)
    EXPECT_DOUBLE_EQ(static_cast<double>(S1[0].real()),
                     static_cast<double>(S2[0].real()));

    // S1(pi) = -S2(pi)
    EXPECT_DOUBLE_EQ(static_cast<double>(S1[1].real()),
                     -static_cast<double>(S2[1].real()));
  }
}
#endif

//******************************************************************************
// TEST(BulkSphere, DISABLED_HandlesInput) {
TEST(BulkSphere, MultiLayerDu) {
  nmie::MultiLayerMie<nmie::FloatType> nmie;
  for (const auto& data : parameters_and_results) {
    auto x = std::get<0>(data);
    auto m = std::get<1>(data);
    //    auto Nstop = nmie::LeRu_near_field_cutoff(m*x)+1;
    nmie.SetLayersSize({x});
    nmie.SetLayersIndex({m});
    //    nmie.SetMaxTerms(Nstop);
    nmie.RunMieCalculation();
    double Qext = static_cast<double>(nmie.GetQext());
    double Qsca = static_cast<double>(nmie.GetQsca());
    double Qext_Du = std::get<2>(data);
    double Qsca_Du = std::get<3>(data);
    EXPECT_FLOAT_EQ(Qext_Du, Qext)
        << "Extinction of the bulk sphere, test case:" << std::get<4>(data)
        << "\nnmax_ = " << nmie.GetMaxTerms();
    EXPECT_FLOAT_EQ(Qsca_Du, Qsca)
        << "Scattering of the bulk sphere, test case:" << std::get<4>(data);
  }
}

//******************************************************************************
TEST(BulkSphere, MesoMieDu) {
  nmie::MultiLayerMie<nmie::FloatType> nmie;
  nmie::MesoMie<nmie::FloatType> mesomie;
  for (const auto& data : parameters_and_results) {
    auto x = std::get<0>(data);
    auto m = std::get<1>(data);
    mesomie.calc_ab(1e-10,    // R
                    x,        // xd
                    x * m,    // xm
                    {1, 0},   // eps_d
                    m * m,    // eps_m
                    {0, 0},   // d_parallel
                    {0, 0});  // d_perp
    mesomie.calc_Q();

    double Qext = static_cast<double>(mesomie.GetQext());
    double Qsca = static_cast<double>(mesomie.GetQsca());
    double Qext_Du = std::get<2>(data);
    double Qsca_Du = std::get<3>(data);
    EXPECT_FLOAT_EQ(Qext_Du, Qext)
        << "Extinction of the bulk sphere, test case:" << std::get<4>(data)
        << "\nnmax_ = " << nmie.GetMaxTerms();
    EXPECT_FLOAT_EQ(Qsca_Du, Qsca)
        << "Scattering of the bulk sphere, test case:" << std::get<4>(data);
  }
}
int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
