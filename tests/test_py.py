import unittest
import numpy as np
from scattnlay import scattnlay, mie, mie_mp

test_cases = [
    # x, {Re(m), Im(m)}, Qext, Qsca, test_name
    [0.099, 0.75+0j, 7.417859e-06, 7.417859e-06],  # ,'a'],
    [0.101, 0.75+0j, 8.033538e-06, 8.033538e-06],  # ,'b'],
    [10, 0.75+0j, 2.232265, 2.232265],  # ,'c'],
    [0.055, 1.5+1j, 0.10149104, 1.131687e-05],  # ,'g'],
    [0.056, 1.5+1j, 0.1033467, 1.216311e-05],  # ,'h'],
    [1, 10+10j, 2.532993, 2.049405],  # ,'k'],
    [100, 1.33+1e-5j, 2.101321, 2.096594],  # ,'e'],
    [100, 1.5+1j, 2.097502, 1.283697],  # ,'i'],
    [1000, 0.75+0j, 1.997908, 1.997908],  # ,'d'],
    [100, 10+10j, 2.071124, 1.836785],  # ,'l'],
    [10000, 1.33+1e-5j, 2.004089, 1.723857],  # ,'f'],
    [10000, 1.5+1j, 2.004368, 1.236574],  # ,'j'],
    [10000, 10+10j, 2.005914, 1.795393],  # ,'m'],
]


class TestStringMethods(unittest.TestCase):

    def test_bulk(self):
        tol = 3e-7
        for solver in [mie, mie_mp]:
            for case in test_cases:
                solver.SetLayersSize(case[0])
                solver.SetLayersIndex(case[1])
                solver.RunMieCalculation()
                Qext = solver.GetQext()
                Qsca = solver.GetQsca()
                self.assertTrue((case[2]-Qext)/Qext < tol)
                self.assertTrue((case[3]-Qsca)/Qsca < tol)


if __name__ == '__main__':
    unittest.main()
