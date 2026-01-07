#!/bin/bash
set -e

# Re-install to ensure fresh build settings
# echo "Ensuring Release build..."
# ./build-all.sh > /dev/null

echo "=========================================================="
echo "          SCATTNLAY PERFORMANCE BENCHMARK (Batch 7)"
echo "=========================================================="
echo ""

# 1. Scalar (Baseline)
# Use smaller N to avoid waiting too long
N_SCALAR=1000
echo "[1] Scalar Baseline (Python Loop + Single Core C++)"
echo "    N=$N_SCALAR"
python3 examples/benchmark-full.py --mode scalar --points $N_SCALAR --repeats 3
echo ""

# 2. SIMD Serial (Batch + Highway, 1 Thread)
N_SIMD=100000
echo "[2] SIMD Serial (Batch + Highway, 1 Thread)"
echo "    N=$N_SIMD"
OMP_NUM_THREADS=1 python3 examples/benchmark-full.py --mode simd --points $N_SIMD --repeats 5
echo ""

# 3. SIMD Parallel (Batch + Highway + OpenMP)
# Use default threads (usually available cores)
echo "[3] SIMD Parallel (Batch + Highway + OpenMP)"
echo "    N=$N_SIMD"
# Unset just in case
unset OMP_NUM_THREADS
# Or force to a high number to verify scaling
export OMP_NUM_THREADS=8
echo "    OMP_NUM_THREADS=$OMP_NUM_THREADS"
python3 examples/benchmark-full.py --mode simd --points $N_SIMD --repeats 10
echo ""

echo "=========================================================="
echo "Benchmark Complete."
