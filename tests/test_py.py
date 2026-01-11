import unittest
from scattnlay import mie, mie_mp, mie_simd
from scattnlay import mesomie
import numpy as np

def round_sig_figs(val, n):
    if val == 0: return 0
    # Format to scientific notation and back to capture precision
    return float('{:g}'.format(float('{:.{p}g}'.format(val, p=n))))

# A list of tests for a bulk sphere from
# Hong Du, "Mie-scattering calculation," Appl. Opt. 43, 1951-1956 (2004)
# table 1: sphere size and refractive index
# followed by resulting extinction and scattering efficiencies
test_cases = [
    # x, m, Qext, Qsca #, test_name
    [0.099, 0.75 + 0j, 7.417859e-06, 7.417859e-06],  # ,'a'],
    [0.101, 0.75 + 0j, 8.033538e-06, 8.033538e-06],  # ,'b'],
    [10, 0.75 + 0j, 2.232265, 2.232265],  # ,'c'],
    [0.055, 1.5 + 1j, 0.10149104, 1.131687e-05],  # ,'g'],
    [0.056, 1.5 + 1j, 0.1033467, 1.216311e-05],  # ,'h'],
    [1, 10 + 10j, 2.532993, 2.049405],  # ,'k'],
    [100, 1.33 + 1e-5j, 2.101321, 2.096594],  # ,'e'],
    [100, 1.5 + 1j, 2.097502, 1.283697],  # ,'i'],
    [1000, 0.75 + 0j, 1.997908, 1.997908],  # ,'d'],
    [100, 10 + 10j, 2.071124, 1.836785],  # ,'l'],
    [10000, 1.33 + 1e-5j, 2.004089, 1.723857],  # ,'f'],
    [10000, 1.5 + 1j, 2.004368, 1.236574],  # ,'j'],
    [10000, 10 + 10j, 2.005914, 1.795393],  # ,'m'],
    # [1.8263846985116234, 0.02867488311561525+1.2957040351341687j, 3, 3]
    # 'water r=1mkm scattnlay 2020/04/22']
    [2.03575204, 1.4558642 + 0.20503704j, 1.952484, 0.9391477],
]


class TestStringMethods(unittest.TestCase):
    def test_bulk_mesomie(self):
        tol = 3e-7
        for solver in [mesomie]:
            print("Using solver: ", solver)
            for case in test_cases:
                x = case[0]
                m = case[1]
                solver.calc_ab(
                    x,  # R
                    x,  # xd
                    x * m,  # xm
                    1,  # eps_d
                    m * m,  # eps_m
                    0,  # d_parallel
                    0,
                )  # d_perp
                solver.calc_Q()
                Qext = solver.GetQext()
                Qsca = solver.GetQsca()
                # print(x, m, Qext)
                self.assertTrue((case[2] - Qext) / Qext < tol)
                self.assertTrue((case[3] - Qsca) / Qsca < tol)

    def test_bulk_multilayer(self):
        tol = 3e-7
        for solver in [mie, mie_mp, mie_simd]:
            if solver is None:
                continue
            print("Using solver: ", solver)
            for case in test_cases:
                solver.SetLayersSize(case[0])
                solver.SetLayersIndex(case[1])
                solver.RunMieCalculation()
                Qext = solver.GetQext()
                Qsca = solver.GetQsca()

                if case == test_cases[0]:
                    print("test case:", case)
                    print("ext tol:", (case[2] - Qext) / Qext)
                    print("sca tol:", (case[3] - Qsca) / Qsca)
                self.assertTrue((case[2] - Qext) / Qext < tol)
                self.assertTrue((case[3] - Qsca) / Qsca < tol)

    def test_du_table_2(self):
        cases = [
            {'x': 0.099, 'm': 0.75 + 0j, 's1_0': 1.81756e-8 + 1.65423e-4j, 's1_pi': 1.81756e-8 + 1.64810e-4j, 'id': 'a'},
            {'x': 0.101, 'm': 0.75 + 0j, 's1_0': 2.04875e-8 + 1.75642e-4j, 's1_pi': 2.04875e-8 + 1.74965e-4j, 'id': 'b'},
            {'x': 10.0, 'm': 0.75 + 0j, 's1_0': 55.8066 + 9.75810j, 's1_pi': -1.07857 + 0.0360881j, 'id': 'c'},
            {'x': 1000.0, 'm': 0.75 + 0j, 's1_0': 499477.0 + 13365.0j, 's1_pi': 17.0578 - 484.251j, 'id': 'd'},
            {'x': 100.0, 'm': 1.33 + 1e-5j, 's1_0': 5253.3 + 124.319j, 's1_pi': -56.5921 - 46.5097j, 'id': 'e'},
            # {'x': 10000.0, 'm': 1.33 + 1e-5j, 's1_0': 5.01022e7 + 153582.0j, 's1_pi': -182.119 + 951.912j, 'id': 'f'},
            {'x': 0.055, 'm': 1.5 + 1j, 's1_0': 7.67526e-5 - 8.34388e-5j, 's1_pi': 7.66140e-5 - 8.33814e-5j, 'id': 'g'},
            {'x': 0.056, 'm': 1.5 + 1j, 's1_0': 8.10238e-5 - 8.80725e-5j, 's1_pi': 8.08721e-5 - 8.80098e-5j, 'id': 'h'},
            {'x': 100.0, 'm': 1.5 + 1j, 's1_0': 5243.75 + 293.417j, 's1_pi': -20.2936 - 4.38444j, 'id': 'i'},
            {'x': 10000.0, 'm': 1.5 + 1j, 's1_0': 5.01092e7 + 175340.0j, 's1_pi': -218.472 + 2064.61j, 'id': 'j'},
            # {'x': 1.0, 'm': 10.0 + 10j, 's1_0': 0.633248 - 0.417931j, 's1_pi': 0.448546 - 0.791236j, 'id': 'k'},
            {'x': 100.0, 'm': 10.0 + 10j, 's1_0': 5177.81 + 26.3381j, 's1_pi': -41.4538 + 18.2181j, 'id': 'l'},
            {'x': 10000.0, 'm': 10.0 + 10j, 's1_0': 5.01479e7 + 120600.0j, 's1_pi': 2252.48 + 3924.47j, 'id': 'm'}
        ]

        for solver in [mie_mp, mie, mie_simd]:
            if solver is None: continue
            print(f"Testing Du Table 2 with solver: {solver}")
            for c in cases:
                if solver == mie_mp and c['x'] > 1000:
                    continue
                solver.SetLayersSize(c['x'])
                solver.SetLayersIndex(c['m'])
                solver.SetAngles(np.array([0.0, np.pi]))
                solver.RunMieCalculation()
                
                S1 = solver.GetS1()
                S2 = solver.GetS2()
                
                # Check S1(0)
                s1_0_calc = S1[0]
                re_0 = round_sig_figs(s1_0_calc.real, 6)
                im_0 = round_sig_figs(s1_0_calc.imag, 6)
                
                self.assertAlmostEqual(re_0, c['s1_0'].real, places=15, msg=f"Case {c['id']} S1(0) real")
                self.assertAlmostEqual(im_0, c['s1_0'].imag, places=15, msg=f"Case {c['id']} S1(0) imag")
                
                # Check S1(pi)
                s1_pi_calc = S1[1]
                re_pi = round_sig_figs(s1_pi_calc.real, 6)
                im_pi = round_sig_figs(s1_pi_calc.imag, 6)
                
                self.assertAlmostEqual(re_pi, c['s1_pi'].real, places=15, msg=f"Case {c['id']} S1(pi) real")
                self.assertAlmostEqual(im_pi, c['s1_pi'].imag, places=15, msg=f"Case {c['id']} S1(pi) imag")

                # S1(0) = S2(0)
                self.assertEqual(S1[0].real, S2[0].real)
                
                # S1(pi) = -S2(pi)
                self.assertEqual(S1[1].real, -S2[1].real)

if __name__ == "__main__":
    unittest.main()
