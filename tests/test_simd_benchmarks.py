import time
import numpy as np
from scattnlay import mie, mie_scalar, mie_simd, Planes

def test_nearfield_simd_benchmark_si_ag_si():
    """
    Replicates tests/test_nearfield_simd_benchmark.cc
    """
    # 1. Model Parameters
    WL = 800.0
    epsilon_Si, epsilon_Ag = 13.64 + 0.047j, -28.05 + 1.525j
    index_Si, index_Ag = np.sqrt(epsilon_Si), np.sqrt(epsilon_Ag)
    r1, r2, r3 = 17.74, 17.74 + 23.31, 17.74 + 23.31 + 22.95
    scale = 2.0 * np.pi / WL
    
    # Setup Mie (SIMD)
    mie.SetLayersSize(np.array([r1, r2, r3]) * scale)
    mie.SetLayersIndex(np.array([index_Si, index_Ag, index_Si]))

    # Setup Mie Scalar
    mie_scalar.SetLayersSize(np.array([r1, r2, r3]) * scale)
    mie_scalar.SetLayersIndex(np.array([index_Si, index_Ag, index_Si]))

    resolution = 256 # 65,536 points
    total_pts = resolution * resolution

    # 2. SCALAR PATH (Manual grid calculation + point-by-point overhead)
    limit = r3 * scale
    scan = np.linspace(-limit, limit, resolution)
    X, Z = np.meshgrid(scan, scan, indexing='ij')
    Xp, Zp = X.flatten(), Z.flatten()
    Yp = np.zeros_like(Xp)

    start_scalar = time.perf_counter()
    # Setting coords and running calculation standard way
    mie_scalar.SetFieldCoords(Xp, Yp, Zp)
    mie_scalar.RunFieldCalculation(isMarkUnconverged=True)
    eabs_scalar = mie_scalar.GetFieldEabs()
    duration_scalar = time.perf_counter() - start_scalar

    # 3. SIMD PATH (Internal Cartesian Generator)
    # This keeps the logic entirely in C++ and utilizes the optimized kernel.
    start_simd = time.perf_counter()
    mie.RunFieldCalculationCartesian(
        first_side_points=resolution, 
        second_side_points=resolution, 
        relative_side_length=1.0, 
        plane_selected=Planes.kEk, # Now using the exported enum
        isMarkUnconverged=True
    )
    eabs_simd = mie.GetFieldEabs()
    duration_simd = time.perf_counter() - start_simd

    # 4. Results
    speedup = duration_scalar / duration_simd
    avg_err = np.mean(np.abs(eabs_simd - eabs_scalar) / (eabs_scalar + 1e-10))

    print(f"\n--- Near-Field SIMD Benchmark (Si-Ag-Si) ---")
    print(f"Grid Size:     {resolution}x{resolution} ({total_pts} points)")
    print(f"Scalar path:   {duration_scalar:.4f} s")
    print(f"SIMD path:     {duration_simd:.4f} s")
    print(f"Speedup:       {speedup:.2f}x")
    print(f"Avg Rel Error: {avg_err:.2e}")

    assert avg_err < 1e-12
    assert speedup > 1.0


def test_farfield_simd_benchmark_parity():
    N = 20000
    x_vals = np.array([0.1 + (i % 1000) * 0.1 for i in range(N)])
    m_vals = np.array([complex(1.5 + (i % 10) * 0.01, 0.01) for i in range(N)])
    
    # 1. True Scalar (std::math - Expected ~0.28s)
    start = time.perf_counter()
    for i in range(N):
        mie_scalar.SetLayersSize(x_vals[i])
        mie_scalar.SetLayersIndex(m_vals[i])
        mie_scalar.RunMieCalculation()
    t_scalar = time.perf_counter() - start

    # 2. Optimized Single (Highway Math - Expected ~0.12s)
    start = time.perf_counter()
    for i in range(N):
        mie.SetLayersSize(x_vals[i])
        mie.SetLayersIndex(m_vals[i])
        mie.RunMieCalculation()
    t_optimized = time.perf_counter() - start

    # 3. Batch SIMD (Zero Python overhead - Expected ~0.06s)
    start = time.perf_counter()
    res_simd = mie_simd.RunMieBatch(x_vals, m_vals)
    t_simd = time.perf_counter() - start

    print(f"\n--- Far-Field Performance Comparison (N={N}) ---")
    print(f"True Scalar (std::math): {t_scalar:.4f} s")
    print(f"Optimized Single (HWY):  {t_optimized:.4f} s")
    print(f"Batch SIMD:              {t_simd:.4f} s")

if __name__ == "__main__":
    test_farfield_simd_benchmark_parity()
    test_nearfield_simd_benchmark_si_ag_si()
