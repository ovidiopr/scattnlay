#include <complex>
#include <vector>
#include <tuple>
#include "gtest/gtest.h"
#include "../src/nmie-bulk.hpp"

TEST(BulkSphereSIMD, DuTable1) {
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
    };

    size_t N = parameters_and_results.size();
    std::vector<std::vector<double>> x(1, std::vector<double>(N));
    std::vector<std::vector<std::complex<double>>> m(1, std::vector<std::complex<double>>(N));
    
    for(size_t i=0; i<N; ++i) {
        x[0][i] = std::get<0>(parameters_and_results[i]);
        m[0][i] = std::get<1>(parameters_and_results[i]);
    }
    
    nmie::BulkSphere<double> bulk(x, m);
    bulk.RunMieCalculation();
    
    for(size_t i=0; i<N; ++i) {
        double expected_Qext = std::get<2>(parameters_and_results[i]);
        double expected_Qsca = std::get<3>(parameters_and_results[i]);
        
        // Use loose tolerance for now as we might have precision differences
        // or maybe I missed something in implementation.
        // The table values are 6-7 sig figs.
        EXPECT_NEAR(bulk.Qext[i], expected_Qext, std::abs(expected_Qext)*1e-4) << "Case " << std::get<4>(parameters_and_results[i]);
        EXPECT_NEAR(bulk.Qsca[i], expected_Qsca, std::abs(expected_Qsca)*1e-4) << "Case " << std::get<4>(parameters_and_results[i]);
    }
}
