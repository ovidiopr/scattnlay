#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <iomanip>
#include "../../src/nmie-basic.hpp"
#include "../../src/nmie-precision.hpp"

using namespace nmie;

void run_benchmark(int nmax, std::complex<double> z, int iterations) {
  std::vector<std::complex<double>> D1(nmax + 1);
  std::vector<std::complex<double>> D3(nmax + 1);
  std::vector<std::complex<double>> Psi(nmax + 1);
  std::vector<std::complex<double>> Zeta(nmax + 1);

  // Warmup
  evalPsiZetaD1D3(z, Psi, Zeta, D1, D3);

  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iterations; ++i) {
    // Prevent optimization
    z += std::complex<double>(1e-15, 1e-15); 
    evalPsiZetaD1D3(z, Psi, Zeta, D1, D3);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;

  double avg_time_us = elapsed.count() / iterations * 1e6;
  // Rough estimate of cycles per element, assuming 3GHz clock. 
  // This is just for relative comparison.
  double cycles_per_element = (elapsed.count() * 3e9) / (iterations * nmax); 

  std::cout << std::setw(10) << nmax 
            << std::setw(25) << z 
            << std::setw(15) << avg_time_us 
            << std::setw(15) << cycles_per_element << "\n";
}

int main() {
  std::cout << "Benchmark: evalPsiZetaD1D3\n";
  std::cout << std::setw(10) << "nmax" 
            << std::setw(25) << "z" 
            << std::setw(15) << "Time (us)" 
            << std::setw(15) << "~Cycles/Elem" << "\n";
  std::cout << std::string(65, '-') << "\n";

  run_benchmark(100, {1.5, 0.1}, 10000);
  run_benchmark(1000, {1.5, 0.1}, 1000);
  run_benchmark(10000, {1.5, 0.1}, 100);
  
  run_benchmark(100, {10.0, 0.1}, 10000);
  run_benchmark(1000, {10.0, 0.1}, 1000);
  
  return 0;
}
