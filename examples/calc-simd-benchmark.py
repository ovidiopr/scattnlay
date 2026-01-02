import numpy as np
import time
try:
    from scattnlay import mie, mie_simd
except ImportError:
    # Fallback for running directly from source/build without installation
    import sys
    import os
    # Try to find the build directory
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../build_simd/scattnlay')))
    try:
        from scattnlay import mie, mie_simd
    except ImportError:
        print("Could not import scattnlay.mie or mie_simd. Ensure the project is built with SIMD support.")
        sys.exit(1)

def run_benchmark():
    if mie_simd is None:
        print("SIMD extension not found. Please build with -DWITH_HWY=ON")
        return

    # --- Configuration ---
    R = 100.0
    m_val = 4.0 + 0.01j
    wl_points = 1000
    wl = np.linspace(400.0, 1000.0, wl_points)
    x = 2.0 * np.pi * R / wl
    m = np.full(wl_points, m_val, dtype=np.complex128)

    print(f"Simulation: {wl_points} points, R={R}nm, m={m_val}")
    print("-----------------------------------------------------")

    # --- 1) Initial Scalar Run (10 times) ---
    def execute_scalar():
        qext, qsca, qabs = [], [], []
        for xi in x:
            mie.SetLayersSize(np.array([xi]))
            mie.SetLayersIndex(np.array([m_val]))
            mie.RunMieCalculation()
            qext.append(mie.GetQext())
            qsca.append(mie.GetQsca())
            qabs.append(mie.GetQabs())
        return np.array(qext), np.array(qsca), np.array(qabs)

    print("Warming up scalar version...")
    t0 = time.perf_counter()
    for _ in range(10):
        res_scalar = execute_scalar()
    t1 = time.perf_counter()
    
    avg_time_per_run = (t1 - t0) / 10.0
    
    # --- 2) Calculate Scaled Repeats for ~5 seconds ---
    total_runs = int(5.0 / avg_time_per_run)
    if total_runs < 1: total_runs = 1
    
    print(f"Running scalar version {total_runs} times to reach ~5s benchmark...")
    
    start_scalar = time.perf_counter()
    for _ in range(total_runs):
        res_scalar = execute_scalar()
    end_scalar = time.perf_counter()
    scalar_duration = end_scalar - start_scalar

    # --- 3) Run SIMD Version (Batch Mode) ---
    print(f"Running SIMD version {total_runs} times...")
    
    start_simd = time.perf_counter()
    for _ in range(total_runs):
        # RunMieBatch accepts the full array of x and m
        res_simd_dict = mie_simd.RunMieBatch(x, m)
    end_simd = time.perf_counter()
    simd_duration = end_simd - start_simd

    # --- 4) Accuracy Comparison ---
    # SIMD output is a dict of lists
    simd_ext = np.array(res_simd_dict["Qext"])
    simd_sca = np.array(res_simd_dict["Qsca"])
    simd_abs = np.array(res_simd_dict["Qabs"])

    def get_diffs(ref, val):
        rel = np.abs(ref - val) / np.where(ref > 0, ref, 1.0)
        return np.mean(rel), np.max(rel)

    avg_e, max_e = get_diffs(res_scalar[0], simd_ext)
    avg_s, max_s = get_diffs(res_scalar[1], simd_sca)
    avg_a, max_a = get_diffs(res_scalar[2], simd_abs)

    # --- 5) Print Results ---
    speedup = scalar_duration / simd_duration

    print("\n" + "="*30)
    print(f"PERFORMANCE")
    print(f"Scalar Duration: {scalar_duration:.4f} s")
    print(f"SIMD Duration:   {simd_duration:.4f} s")
    print(f"Speedup Factor:  {speedup:.2f}x")
    print("="*30)
    print(f"ACCURACY (Relative Difference)")
    print(f"Qext: Avg={avg_e:.2e}, Max={max_e:.2e}")
    print(f"Qsca: Avg={avg_s:.2e}, Max={max_s:.2e}")
    print(f"Qabs: Avg={avg_a:.2e}, Max={max_a:.2e}")
    print("="*30)

if __name__ == "__main__":
    run_benchmark()
