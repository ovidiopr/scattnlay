#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>
#include "../src/nmie-basic.hpp"
#include "../src/nmie-nearfield.hpp"

using namespace nmie;

// Result container for the benchmark
struct NearFieldBenchResult {
  std::chrono::duration<double> duration;
  std::vector<double> Eabs;
};

/**
 * Template helper to execute the Near-Field Cartesian calculation for a
 * multilayer sphere.
 */
template <typename Engine>
NearFieldBenchResult measureNearFieldCartesian(
    const std::vector<double>& x,
    const std::vector<std::complex<double>>& m,
    int res) {
  MultiLayerMie<double, Engine> nmie_obj;
  nmie_obj.SetLayersSize(x);
  nmie_obj.SetLayersIndex(m);

  auto start = std::chrono::high_resolution_clock::now();

  // 512x512 grid, Ek plane, relative side length of 2.0 (covers 2x total
  // radius)
  nmie_obj.RunFieldCalculationCartesian(res, res, 2.0, Planes::kEk, 0.0, 0.0,
                                        0.0, true);

  std::vector<double> results = nmie_obj.GetFieldEabs();

  auto end = std::chrono::high_resolution_clock::now();
  return {end - start, std::move(results)};
}

TEST(SIMDBenchmark, NearFieldMultilayerSiAgSi) {
#ifndef WITH_HWY
  GTEST_SKIP() << "Skipping benchmark: SIMD (Highway) not enabled.";
#endif

  // Model parameters provided
  const double WL = 800.0;  // nm
  const std::complex<double> epsilon_Si = {13.64, 0.047};
  const std::complex<double> epsilon_Ag = {-28.05, 1.525};
  const std::complex<double> index_Si = std::sqrt(epsilon_Si);
  const std::complex<double> index_Ag = std::sqrt(epsilon_Ag);

  const double core_width = 17.74;   // nm
  const double inner_width = 23.31;  // nm
  const double outer_width = 22.95;  // nm

  // Cumulative radii for multilayer size parameters
  double r1 = core_width;
  double r2 = r1 + inner_width;
  double r3 = r2 + outer_width;

  // Convert to size parameters: x = 2 * pi * r / WL
  const double scale = 2.0 * std::numbers::pi / WL;
  std::vector<double> x_vals = {r1 * scale, r2 * scale, r3 * scale};
  std::vector<std::complex<double>> m_vals = {index_Si, index_Ag, index_Si};

  const int resolution = 256;

  // 1. Run Benchmark via Helper for both Engines
  auto scalar = measureNearFieldCartesian<ScalarEngine<double>>(x_vals, m_vals,
                                                                resolution);
  auto simd = measureNearFieldCartesian<DefaultEngine<double>>(x_vals, m_vals,
                                                               resolution);

  // 2. Calculate Parity (ignoring NaNs to ensure robustness)
  double sum_rel_err = 0;
  size_t valid_pts = 0;
  size_t total_pts = scalar.Eabs.size();

  for (size_t i = 0; i < total_pts; ++i) {
    double s_val = scalar.Eabs[i];
    double v_val = simd.Eabs[i];

    if (std::isfinite(s_val) && std::isfinite(v_val) && s_val > 1e-10) {
      sum_rel_err += std::abs(v_val - s_val) / s_val;
      valid_pts++;
    }
  }

  double avg_err = (valid_pts > 0) ? (sum_rel_err / valid_pts) : 0.0;
  double speedup = scalar.duration.count() / simd.duration.count();

  // 3. Print Results
  std::cout << "--- Near-Field Multilayer Benchmark (Si-Ag-Si) ---"
            << std::endl;
  std::cout << "Grid Size:     " << resolution << "x" << resolution
            << std::endl;
  std::cout << "Scalar time:   " << std::fixed << std::setprecision(4)
            << scalar.duration.count() << " s" << std::endl;
  std::cout << "SIMD time:     " << simd.duration.count() << " s" << std::endl;
  std::cout << "Speedup:       " << std::setprecision(2) << speedup << "x"
            << std::endl;
  std::cout << "Avg Rel Error: " << std::scientific << std::setprecision(2)
            << avg_err << " (" << valid_pts << " valid points)" << std::endl;

  // 4. Assertions
  EXPECT_LT(avg_err, 1e-12);
  EXPECT_GT(speedup, 1.0);
  EXPECT_EQ(valid_pts, total_pts);  // Ensure high convergence rate
}