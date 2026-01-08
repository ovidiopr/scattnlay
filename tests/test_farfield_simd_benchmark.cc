#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "../src/nmie-basic.hpp"
#include "../src/nmie-batch.hpp"

using namespace nmie;

struct FarFieldResult {
    std::chrono::duration<double> duration;
    std::vector<double> Qext;
};

// Template helper to measure Far-Field performance
template <typename Engine>
FarFieldResult measureFarField(const MieBatchInput& input) {
    auto start = std::chrono::high_resolution_clock::now();
    size_t N = input.x.size();
    std::vector<double> qext_results(N);

    if constexpr (std::is_same_v<Engine, ScalarEngine<double>>) {
      MultiLayerMie<double, ScalarEngine<double>> scalar_mie;
      for (size_t i = 0; i < N; ++i) {
        scalar_mie.SetLayersSize({input.x[i]});
        scalar_mie.SetLayersIndex({input.m[i]});
        scalar_mie.RunMieCalculation();
        qext_results[i] = scalar_mie.GetQext();
      }
    } else {
      auto output = RunMieBatch<double>(input);
      qext_results = std::move(output.Qext);
    }

    auto end = std::chrono::high_resolution_clock::now();
    return {end - start, std::move(qext_results)};
}

TEST(SIMDBenchmark, FarFieldParityAndSpeedup) {
#ifndef WITH_HWY
    GTEST_SKIP() << "Skipping benchmark: SIMD (Highway) not enabled.";
#endif

    const int N = 20000;
    MieBatchInput input;
    input.x.reserve(N);
    input.m.reserve(N);

    for (int i = 0; i < N; ++i) {
        input.x.push_back(0.1 + (i % 1000) * 0.1);
        input.m.push_back({1.5 + (i % 10) * 0.01, 0.01});
    }

    // 1. Scalar Reference (1-core)
    auto scalar = measureFarField<ScalarEngine<double>>(input);

    // Helper for reporting per-core results
    auto report = [&](const std::string& label, const FarFieldResult& res) {
      double sum_rel_err = 0;
      for (int i = 0; i < N; ++i) {
        sum_rel_err += std::abs(res.Qext[i] - scalar.Qext[i]) / scalar.Qext[i];
      }
      std::cout << label << ":" << std::endl;
      std::cout << "  - time:   " << res.duration.count() << " s" << std::endl;
      std::cout << "  - speedup: "
                << scalar.duration.count() / res.duration.count() << "x"
                << std::endl;
      std::cout << "  - avg error:   " << sum_rel_err / N << std::endl;
    };

    // 2. SIMD 1-core
    nmie::setNumThreads(1);
    auto simd1 = measureFarField<DefaultEngine<double>>(input);

    // 3. SIMD Multi-core
    int num_cores = 12;
    nmie::setNumThreads(num_cores);
    auto simdn = measureFarField<DefaultEngine<double>>(input);

    // Output formatting
    std::cout << "--- Far-Field SIMD Benchmark (N=" << N << ") ---" << std::endl;
    std::cout << "Scalar time: " << scalar.duration.count() << " s" << std::endl;

    report("SIMD 1-core", simd1);
    report("SIMD " + std::to_string(num_cores) + "-core", simdn);

    // Validation
    EXPECT_LT(std::abs(simd1.Qext[0] - scalar.Qext[0]), 1e-12);
}