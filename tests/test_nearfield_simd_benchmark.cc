#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
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
 * Template helper to execute the Near-Field Cartesian calculation.
 * Replicates the logic of the Python benchmark script:
 * 512x512 grid, Ek plane, relative side length of 2.0.
 */
template <typename Engine>
NearFieldBenchResult measureNearFieldCartesian(double x,
                                               std::complex<double> m,
                                               int res) {
  MultiLayerMie<double, Engine> nmie_obj;
  nmie_obj.SetLayersSize({x});
  nmie_obj.SetLayersIndex({m});

  auto start = std::chrono::high_resolution_clock::now();

  // Exactly replicating the Python benchmark call:
  // resolution, resolution, side_length, plane, at_x, at_y, at_z,
  // markUnconverged
  nmie_obj.RunFieldCalculationCartesian(res, res, 2.0, Planes::kEk, 0.0, 0.0,
                                        0.0, true);

  std::vector<double> results = nmie_obj.GetFieldEabs();

  auto end = std::chrono::high_resolution_clock::now();
  return {end - start, std::move(results)};
}

TEST(SIMDBenchmark, NearFieldParityAndSpeedup) {
#ifndef WITH_HWY
  GTEST_SKIP() << "Skipping benchmark: SIMD (Highway) not enabled.";
#endif

  // Setup identical to examples/calc-simd-benchmark-nearfield.py
  const double x = 20.0;
  const std::complex<double> m = {1.5, 0.1};
  const int resolution = 512;

  // 1. Run Benchmark via Helper
  auto scalar =
      measureNearFieldCartesian<ScalarEngine<double>>(x, m, resolution);
  auto simd =
      measureNearFieldCartesian<DefaultEngine<double>>(x, m, resolution);

  // 2. Calculate Parity (ignoring NaNs to match benchmark behavior)
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
  std::cout << "--- Near-Field SIMD Benchmark (Cartesian " << resolution << "x"
            << resolution << ") ---" << std::endl;
  std::cout << "Scalar time:   " << std::fixed << std::setprecision(4)
            << scalar.duration.count() << " s" << std::endl;
  std::cout << "SIMD time:     " << simd.duration.count() << " s" << std::endl;
  std::cout << "Speedup:       " << std::setprecision(2) << speedup << "x"
            << std::endl;
  std::cout << "Avg Rel Error: " << std::scientific << std::setprecision(2)
            << avg_err << " (" << valid_pts << " valid points)" << std::endl;

  // 4. Assertions
  // We expect extremely high parity (1e-12 or better)
  EXPECT_LT(avg_err, 1e-12);

  // We expect at least some speedup on any SIMD hardware
  EXPECT_GT(speedup, 1.0);

  // Verify that the grid is mostly convergent
  EXPECT_GT(valid_pts, total_pts * 0.999);
}