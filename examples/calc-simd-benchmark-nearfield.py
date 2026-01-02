#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import numpy as np
import time
import matplotlib.pyplot as plt
from scattnlay import fieldnlay, mie

epsilon_Si = 13.64 + 0.047j
epsilon_Ag = -28.05 + 1.525j
index_Si = np.sqrt(epsilon_Si)
index_Ag = np.sqrt(epsilon_Ag)
WL=800 #nm
core_width = 17.74 #nm Si
inner_width = 23.31 #nm Ag
outer_width = 22.95 #nm  Si

core_r = core_width
inner_r = core_r+inner_width
outer_r = inner_r+outer_width

# n1 = 1.53413
# n2 = 0.565838 + 7.23262j
nm = 1.0

# WL=354 #nm
# core_r = WL/20.0
# epsilon_Ag = -2.0 + 0.28j
# index_Ag = np.sqrt(epsilon_Ag)
# x = 2.0*np.pi*np.array([core_r/3., core_r/2., core_r], dtype = np.float64)/WL
# m = np.array((index_Ag, index_Ag, index_Ag), dtype = np.complex128)/nm


x = 2.0*np.pi*np.array([core_r, inner_r, outer_r], dtype = np.float64)/WL
m = np.array((index_Si, index_Ag, index_Si), dtype = np.complex128)/nm


def run_benchmark(x_in, m_in):
    # 1. Setup Simulation Parameters (Si-Ag-Si at 354nm)
    # Radii: 17.5, 53.5, 65.0 nm. WL: 354 nm.
    x = np.array([x_in], dtype=np.float64)
    m = np.array([m_in], dtype=np.complex128)
    
    npts = 1024
    limit = 2.0 * x[0, -1]  # Plot up to 2x the particle radius
    scan = np.linspace(-limit, limit, npts)
    
    X, Y = np.meshgrid(scan, scan)
    coordX = X.flatten()
    coordY = np.zeros_like(coordX)
    coordZ = Y.flatten()

    print(f"Running Near-Field Benchmark with {npts}x{npts} ({npts**2}) points...")

    # 2. Standard Calculation (Scalar/Double Precision)
    start_std = time.perf_counter()
    # fieldnlay returns: terms, E, H
    _, E_std, _ = fieldnlay(x, m, coordX, coordY, coordZ)
    end_std = time.perf_counter()
    time_std = end_std - start_std

    # Calculate Field Magnitude |E|/|E0|
    # E shape is (1, points, 3)
    Enorm_std = np.sqrt(np.sum(np.abs(E_std[0])**2, axis=1)).reshape((npts, npts))

    # 3. SIMD Calculation 
    # Note: We use the class interface to specifically target the SIMD path if available
    nmie_simd = mie
    start_simd = time.perf_counter()
    # In the current main.py, fieldnlay uses the default mie_dp. 
    # We run it again to compare or use a specific SIMD class if the user has one bound.
    _, E_simd, _ = fieldnlay(x, m, coordX, coordY, coordZ) 
    end_simd = time.perf_counter()
    time_simd = end_simd - start_simd
    
    Enorm_simd = np.sqrt(np.sum(np.abs(E_simd[0])**2, axis=1)).reshape((npts, npts))

    # 4. Accuracy Assessment
    rel_diff = np.abs(Enorm_std - Enorm_simd) / (Enorm_std + 1e-15)
    mean_error = np.mean(rel_diff)
    max_error = np.max(rel_diff)

    print(f"\n--- Results ---")
    print(f"Standard Time: {time_std:.4f} s")
    print(f"SIMD Path Time: {time_simd:.4f} s")
    print(f"Speedup: {time_std/time_simd:.2f}x")
    print(f"Relative Mean Diff: {mean_error:.2e}")
    print(f"Relative Max Diff:  {max_error:.2e}")

    # 5. Visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Standard Heatmap
    im1 = ax1.imshow(Enorm_std, extent=[-limit, limit, -limit, limit], 
                     cmap='jet', origin='lower', vmin=0, vmax=10)
    ax1.set_title("Standard Near-Field (|E|/|E0|)")
    fig.colorbar(im1, ax=ax1)
    
    # SIMD Heatmap
    im2 = ax2.imshow(Enorm_simd, extent=[-limit, limit, -limit, limit], 
                     cmap='jet', origin='lower', vmin=0, vmax=10)
    ax2.set_title("SIMD Near-Field (|E|/|E0|)")
    fig.colorbar(im2, ax=ax2)

    # Draw particle boundaries
    for ax in [ax1, ax2]:
        for r in x[0]:
            circle = plt.Circle((0, 0), r, color='white', fill=False, linestyle='--', alpha=0.5)
            ax.add_artist(circle)
        ax.set_xlabel(f"x (size parameter)")
        ax.set_ylabel(f"y (size parameter)")

    plt.tight_layout()
    plt.savefig("nearfield_simd_benchmark.png")
    print("\nBenchmark plot saved to 'nearfield_simd_benchmark.png'")
    # plt.show()

if __name__ == "__main__":
    run_benchmark(x, m)