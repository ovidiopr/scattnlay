#include <complex>
#include <vector>
#include <tuple>
#include "../src/nmie-basic.hpp"
#include "gtest/gtest.h"

TEST(StatelessCoeffs, DuTable1) {
  // x, {Re(m), Im(m)}, Qext, Qsca, test_name
  std::vector<std::tuple<double, std::complex<double>, double, double, char>>
      parameters_and_results{
          {0.099, {0.75, 0}, 7.417859e-06, 7.417859e-06, 'a'},
          {0.101, {0.75, 0}, 8.033538e-06, 8.033538e-06, 'b'},
          {10, {0.75, 0}, 2.232265, 2.232265, 'c'},
          // {1000, {0.75, 0}, 1.997908, 1.997908, 'd'}, // Skip large ones for speed
          {100, {1.33, 1e-5}, 2.101321, 2.096594, 'e'},
          // {10000, {1.33, 1e-5}, 2.004089, 1.723857, 'f'},
          {0.055, {1.5, 1}, 0.101491, 1.131687e-05, 'g'},
          {0.056, {1.5, 1}, 0.103347, 1.216311e-05, 'h'},
          {100, {1.5, 1}, 2.097502, 1.283697, 'i'},
          // {10000, {1.5, 1}, 2.004368, 1.236574, 'j'},
          {1, {10, 10}, 2.532993, 2.049405, 'k'},
          {100, {10, 10}, 2.071124, 1.836785, 'l'},
          // {10000, {10, 10}, 2.005968, 1.795393, 'm'},
      };

  for (const auto& params : parameters_and_results) {
    double x_val = std::get<0>(params);
    std::complex<double> m_val = std::get<1>(params);

    nmie::MultiLayerMie<double> nmie;
    nmie.SetLayersSize({x_val});
    nmie.SetLayersIndex({m_val});
    nmie.RunMieCalculation();
    
    auto an_ref = nmie.GetAn();
    auto bn_ref = nmie.GetBn();
    int nmax = nmie.GetMaxTerms();

    // Stateless calculation
    nmie::MieBuffers<double, nmie::ScalarEngine<double>> buffers;
    buffers.resize(nmax, 1);
    buffers.updateSize(nmax, 1);

    auto get_x = [&](int l) { return x_val; };
    auto get_m = [&](int l) { return m_val; };

    nmie::calcScattCoeffsKernel<double, nmie::ScalarEngine<double>>(
        nmax, 1, -1, get_x, get_m, buffers);

    for (int n = 0; n < nmax; ++n) {
      EXPECT_NEAR(buffers.an[n].real(), an_ref[n].real(), 1e-14);
      EXPECT_NEAR(buffers.an[n].imag(), an_ref[n].imag(), 1e-14);
      EXPECT_NEAR(buffers.bn[n].real(), bn_ref[n].real(), 1e-14);
      EXPECT_NEAR(buffers.bn[n].imag(), bn_ref[n].imag(), 1e-14);
    }
  }
}
