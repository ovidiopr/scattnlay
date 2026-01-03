#include <complex>
#include <vector>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "../src/nmie-basic.hpp"
#include "../src/nmie-nearfield.hpp"
#include "gtest/gtest.h"

TEST(OpenMPField, ParallelExecution) {
  double x = 10.0;
  std::complex<double> m(1.33, 0.0);
  int nmax = 100; // Fixed nmax for simplicity

  // Setup inputs
  std::vector<double> size_param = {x};
  std::vector<std::complex<double>> refractive_index = {m};
  
  // 1. Calculate Scattering Coeffs (an, bn)
  nmie::MieBuffers<double, nmie::ScalarEngine<double>> master_buffers;
  master_buffers.resize(nmax, 1, 0);
  master_buffers.updateSize(nmax, 1, 0);

  auto get_x = [&](int l) { return x; };
  auto get_m = [&](int l) { return m; };

  nmie::calcScattCoeffsKernel<double, nmie::ScalarEngine<double>>(
      nmax, 1, -1, get_x, get_m, master_buffers);

  // 2. Calculate Expansion Coeffs (aln, bln...)
  nmie::calcExpanCoeffsKernel<double, nmie::ScalarEngine<double>>(
      nmax, refractive_index, size_param, -1, master_buffers);

  // 3. Prepare points
  int num_points = 10000;
  std::vector<double> Xp(num_points), Yp(num_points), Zp(num_points);
  for(int i=0; i<num_points; ++i) {
      Xp[i] = 1.5 * x; // Outside
      Yp[i] = 0.0;
      Zp[i] = 0.0;
  }

  // Outputs
  std::vector<std::vector<std::complex<double>>> Es(num_points), Hs(num_points);
  std::vector<std::vector<double>> coords_polar(num_points);
  for(int i=0; i<num_points; ++i) {
      Es[i].resize(3);
      Hs[i].resize(3);
      coords_polar[i].resize(3);
  }

  // 4. Parallel Field Calculation
  #pragma omp parallel
  {
      // Thread-local buffers
      nmie::MieBuffers<double, nmie::ScalarEngine<double>> thread_buffers;
      thread_buffers.resize(nmax, 1, 0); // Resize for scratch space
      
      #pragma omp for
      for (int i = 0; i < num_points; ++i) {
          // Dummy loop to ensure syntax correctness of pragma
      }
  }
  
  // Better approach: Chunking
#ifdef _OPENMP
  int num_threads = omp_get_max_threads();
#else
  int num_threads = 1;
#endif
  int chunk_size = (num_points + num_threads - 1) / num_threads;

  #pragma omp parallel
  {
      nmie::MieBuffers<double, nmie::ScalarEngine<double>> thread_buffers;
      thread_buffers.resize(nmax, 1, 0);

      #pragma omp for schedule(static)
      for (int i = 0; i < num_points; i += 100) { // Process in blocks of 100
          int count = std::min(100, num_points - i);
          
          nmie::RunFieldKernel<double, nmie::ScalarEngine<double>>(
              Xp, Yp, Zp,
              i, count,
              nmax,
              size_param,
              refractive_index,
              master_buffers.aln,
              master_buffers.bln,
              master_buffers.cln,
              master_buffers.dln,
              thread_buffers,
              Es, Hs,
              coords_polar
          );
      }
  }

  // Verify results (check first point)
  EXPECT_NE(std::abs(Es[0][0]) + std::abs(Es[0][1]) + std::abs(Es[0][2]), 0.0);
}
