#!/usr/bin/env python3
import argparse
import numpy as np
import time
import os
import sys

# Try to import scattnlay package
try:
    from scattnlay import mie as mie_scalar
    from scattnlay import mie_simd
except ImportError:
    # If running from source/build directory
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../build/scattnlay')))
    from scattnlay import mie as mie_scalar
    from scattnlay import mie_simd

def run_scalar(x, m, iterations=1):
    """Run the scalar loop implementation."""
    t0 = time.perf_counter()
    
    for _ in range(iterations):
        qext, qsca, qabs = [], [], []
        for i in range(len(x)):
            mie_scalar.SetLayersSize(np.array([x[i]]))
            mie_scalar.SetLayersIndex(np.array([m[i]]))
            mie_scalar.RunMieCalculation()
            qext.append(mie_scalar.GetQext())
            qsca.append(mie_scalar.GetQsca())
            qabs.append(mie_scalar.GetQabs())
            
    t1 = time.perf_counter()
    return t1 - t0

def run_simd(x, m, theta=None, iterations=1):
    """Run the SIMD batch implementation."""
    # Ensure contiguous C-style arrays
    x_c = np.ascontiguousarray(x, dtype=np.float64)
    m_c = np.ascontiguousarray(m, dtype=np.complex128)
    theta_c = np.ascontiguousarray(theta, dtype=np.float64) if theta is not None else np.array([], dtype=np.float64)

    t0 = time.perf_counter()
    for _ in range(iterations):
        res = mie_simd.RunMieBatch(x_c, m_c, theta_c)
    t1 = time.perf_counter()
    return t1 - t0

def main():
    parser = argparse.ArgumentParser(description="Scattnlay Benchmark Tool")
    parser.add_argument("--mode", choices=["scalar", "simd"], required=True, help="Execution mode")
    parser.add_argument("--points", type=int, default=100000, help="Number of particles")
    parser.add_argument("--repeats", type=int, default=5, help="Measurement repeats")
    
    args = parser.parse_args()

    # --- Setup Problem ---
    # Generate random particle configurations
    np.random.seed(42)
    x = np.random.uniform(0.1, 100.0, args.points)
    m = np.random.uniform(1.33, 1.8, args.points) + 1j * np.random.uniform(0.001, 0.1, args.points)
    
    print(f"Problem Size: {args.points} particles")
    print(f"Mode: {args.mode}")
    if args.mode == "simd":
        print(f"Threads: {os.environ.get('OMP_NUM_THREADS', 'Default')}")

    # --- Warmup ---
    # Run a small batch to warm up caches/JIT
    warmup_x = x[:100]
    warmup_m = m[:100]
    if args.mode == "scalar":
        run_scalar(warmup_x, warmup_m, 1)
    else:
        run_simd(warmup_x, warmup_m, None, 1)

    # --- Benchmark ---
    durations = []
    for i in range(args.repeats):
        # We process the FULL array 'args.points' in one go (or one loop for scalar)
        # This measures the time to process 'points' particles.
        dt = 0
        if args.mode == "scalar":
            dt = run_scalar(x, m, 1)
        else:
            dt = run_simd(x, m, None, 1)
        durations.append(dt)
        print(f"Run {i+1}: {dt:.4f} s")

    avg_time = np.mean(durations)
    std_time = np.std(durations)
    
    print("-" * 40)
    print(f"Average Time: {avg_time:.4f} ± {std_time:.4f} s")
    print(f"Rate: {args.points / avg_time:.2e} calculations/sec")
    print("-" * 40)

if __name__ == "__main__":
    main()
