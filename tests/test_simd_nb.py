try:
    from scattnlay import scattnlay_simd as scs
except ImportError:
    import scattnlay_simd as scs
import numpy as np
import unittest

class TestSimdBatch(unittest.TestCase):
    def test_simd_batch(self):
        # 1. Create large batch
        N = 1000
        x = np.linspace(1.0, 100.0, N, dtype=np.float64)
        m = np.full(N, 1.5 + 0.1j, dtype=np.complex128)
        theta = np.linspace(0, np.pi, 180, dtype=np.float64)

        # 2. Run Batch
        results = scs.RunMieBatch(x, m, theta)

        # 3. Verify types and zero-copy
        for key in ["Qext", "Qsca", "S1"]:
            arr = results[key]
            self.assertIsInstance(arr, np.ndarray, f"{key} is not an ndarray")
            self.assertIsNotNone(arr.base, f"{key} did not use zero-copy (capsule missing)")
        
        # 4. Verify Shape of S1 (N spheres x 180 angles)
        self.assertEqual(results["S1"].shape, (1000, 180))
        print("SIMD Nanobind Migration: SUCCESS")

if __name__ == "__main__":
    unittest.main()
