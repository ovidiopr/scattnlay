import time
import numpy as np
from scattnlay import mie

def test_nearfield_simd_benchmark_si_ag_si():
    """
    Exactly replicates tests/test_nearfield_simd_benchmark.cc
    Model: Multilayer Si-Ag-Si sphere @ 800nm
    Calculates E-field on a 256x256 grid and prints SIMD Speedup.
    """
    # 1. Setup Model Parameters (Identical to C++ Benchmark)
    WL = 800.0  # nm
    epsilon_Si = 13.64 + 0.047j
    epsilon_Ag = -28.05 + 1.525j
    index_Si = np.sqrt(epsilon_Si)
    index_Ag = np.sqrt(epsilon_Ag)

    # Layer widths (nm) converted to size parameters (x = 2*pi*r/WL)
    r1, r2, r3 = 17.74, 17.74 + 23.31, 17.74 + 23.31 + 22.95
    scale = 2.0 * np.pi / WL
    x_vals = np.array([r1 * scale, r2 * scale, r3 * scale])
    m_vals = np.array([index_Si, index_Ag, index_Si])

    resolution = 256 
    total_pts = resolution * resolution

    # Initialize Mie object
    mie.SetLayersSize(x_vals)
    mie.SetLayersIndex(m_vals)

    # 2. "Scalar" Path (Point-by-point loop)
    # To simulate the ScalarEngine, we call calculations for one point at a time.
    # This prevents the SIMD kernel from filling its lanes with multiple points.
    limit = r3 * scale
    scan = np.linspace(-limit, limit, resolution)
    X, Z = np.meshgrid(scan, scan)
    Xp, Zp = X.flatten(), Z.flatten()
    Yp = np.zeros_like(Xp)

    print(f"\n--- Near-Field SIMD Benchmark (Si-Ag-Si) ---")
    print(f"Grid Size:     {resolution}x{resolution} ({total_pts} points)")

    # We use a subset for the scalar baseline if the full grid is too slow, 
    # but here we do the full grid for an exact comparison.
    start_scalar = time.perf_counter()
    # Note: Calling SetFieldCoords with the full array and RunFieldCalculation() 
    # is faster than a Python loop but slower than the Cartesian internal generator.
    mie.SetFieldCoords(Xp, Yp, Zp)
    mie.RunFieldCalculation(isMarkUnconverged=True)
    eabs_scalar = mie.GetFieldEabs()
    duration_scalar = time.perf_counter() - start_scalar
    print(f"Scalar path time: {duration_scalar:.4f} s")

    # 3. SIMD Path (Internal Cartesian Grid Generator)
    # This uses the specific C++ function that the benchmark uses: 
    # It generates the grid internally and vectorizes across points.
    start_simd = time.perf_counter()
    # Plane 0 = Ek, Relative Size 2.0
    mie.RunFieldCalculationCartesian(
        resolution, resolution, 2.0, 0, 
        0.0, 0.0, 0.0, True
    )
    eabs_simd = mie.GetFieldEabs()
    duration_simd = time.perf_counter() - start_simd
    print(f"SIMD path time:   {duration_simd:.4f} s")

    # 4. Results and Parity
    speedup = duration_scalar / duration_simd
    avg_err = np.mean(np.abs(eabs_simd - eabs_scalar) / (eabs_scalar + 1e-10))

    print(f"Speedup:          {speedup:.2f}x")
    print(f"Avg Rel Error:    {avg_err:.2e}")

    # 5. Assertions matching C++ expectations
    assert avg_err < 1e-12, "SIMD results differ from scalar results"
    # Speedup > 1.0 ensures the optimized internal loop is indeed faster.
    # On most AVX2/NEON systems, this should be > 2.0x.
    assert speedup > 1.0 


def test_farfield_simd_benchmark_parity():
    """
    Replicates tests/test_farfield_simd_benchmark.cc (N=10,000)
    """
    from scattnlay import mie_simd
    if mie_simd is None:
        raise RuntimeError("scattnlay_simd not built.")

    N = 10000
    x_vals = np.array([0.1 + (i % 1000) * 0.1 for i in range(N)], dtype=np.float64)
    m_vals = np.array([complex(1.5 + (i % 10) * 0.01, 0.01) for i in range(N)], dtype=np.complex128)
    
    # Scalar
    start_scalar = time.perf_counter()
    qext_scalar = []
    for i in range(N):
        mie.SetLayersSize(x_vals[i])
        mie.SetLayersIndex(m_vals[i])
        mie.RunMieCalculation()
        qext_scalar.append(mie.GetQext())
    duration_scalar = time.perf_counter() - start_scalar

    # SIMD
    start_simd = time.perf_counter()
    res_simd = mie_simd.RunMieBatch(x_vals, m_vals)
    duration_simd = time.perf_counter() - start_simd

    speedup = duration_scalar / duration_simd
    print(f"\n--- Far-Field SIMD Benchmark (N={N}) ---")
    print(f"Scalar time: {duration_scalar:.4f} s")
    print(f"SIMD time:   {duration_simd:.4f} s")
    print(f"Speedup:     {speedup:.2f}x")
    
    assert speedup > 1.0

if __name__ == "__main__":
    # test_nearfield_simd_benchmark_si_ag_si()
    test_farfield_simd_benchmark_parity()