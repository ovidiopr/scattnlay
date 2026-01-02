#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include <chrono>
#include <iostream>
#include "../src/nmie-basic.hpp"
#include "../src/nmie-batch.hpp"

using namespace nmie;

struct FarFieldResult {
    std::chrono::duration<double> duration;
    std::vector<double> Qext;
};

// Template helper to measure Far-Field performance
// For ScalarEngine, it performs a loop. For HighwayEngine, it uses the optimized batch runner.
template <typename Engine>
FarFieldResult measureFarField(const MieBatchInput& input) {
    auto start = std::chrono::high_resolution_clock::now();
    size_t N = input.x.size();
    std::vector<double> qext_results(N);

    if constexpr (std::is_same_v<Engine, ScalarEngine<double>>) {
        // Reference implementation: Loop over individual particles
        MultiLayerMie<double, ScalarEngine<double>> scalar_mie;
        for (size_t i = 0; i < N; ++i) {
            scalar_mie.SetLayersSize({input.x[i]});
            scalar_mie.SetLayersIndex({input.m[i]});
            scalar_mie.RunMieCalculation();
            qext_results[i] = scalar_mie.GetQext();
        }
    } else {
        // Optimized implementation: Batch processing via Google Highway
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

    const int N = 10000;
    MieBatchInput input;
    input.x.reserve(N);
    input.m.reserve(N);

    for (int i = 0; i < N; ++i) {
        input.x.push_back(0.1 + (i % 1000) * 0.1);
        input.m.push_back({1.5 + (i % 10) * 0.01, 0.01});
    }

    auto scalar = measureFarField<ScalarEngine<double>>(input);
    auto simd = measureFarField<DefaultEngine<double>>(input);

    double sum_rel_err = 0;
    for (int i = 0; i < N; ++i) {
        sum_rel_err += std::abs(simd.Qext[i] - scalar.Qext[i]) / scalar.Qext[i];
    }
    double avg_err = sum_rel_err / N;
    double speedup = scalar.duration.count() / simd.duration.count();

    std::cout << "--- Far-Field SIMD Benchmark (N=" << N << ") ---" << std::endl;
    std::cout << "Scalar time: " << scalar.duration.count() << " s" << std::endl;
    std::cout << "SIMD time:   " << simd.duration.count() << " s" << std::endl;
    std::cout << "Speedup:     " << speedup << "x" << std::endl;
    std::cout << "Avg Error:   " << avg_err << std::endl;

    EXPECT_LT(avg_err, 1e-13);
    EXPECT_GT(speedup, 0.9);
}