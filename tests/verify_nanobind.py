from scattnlay import scattnlay_dp as sc
import numpy as np
import sys

def verify_zero_copy():
    mie = sc.mie_dp()
    mie.SetLayersSize(np.array([100.0]))
    mie.SetLayersIndex(np.array([1.5 + 0.1j]))
    mie.RunMieCalculation()

    s1 = mie.GetS1()
    # Check 1: Is it a NumPy array?
    assert isinstance(s1, np.ndarray), "Failed: Not a numpy array"
    
    # Check 2: Does it own its memory? 
    # If zero-copy worked, s1.base should be an nb::capsule (not None)
    assert s1.base is not None, "Failed: Array is a copy, not a view of C++ memory"
    
    print("Zero-copy verification: PASSED")

if __name__ == "__main__":
    verify_zero_copy()
