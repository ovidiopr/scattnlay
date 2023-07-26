from scattnlay import mesomie, mie
import numpy as np



def eval_mie_single_spectrum_point(x,m,
                                   R=None, d_perp=None, d_parl=None):
    mie.SetLayersSize(x)
    mie.SetLayersIndex(m)
    mie.RunMieCalculation()
    return mie.GetQext(), mie.GetQsca()


def eval_mesomie_single_spectrum_point(x, m, n_m, R, d_perp, d_parl):
    mesomie.calc_ab(R*10,      # calc_ab needs R in angstrem
                x,      # xd
                x * m,  # xm
                n_m * n_m,      # eps_d 
                m * m,  # eps_m
                d_parl,      # d_parallel
                d_perp)      # d_perp
    mesomie.calc_Q()
    return mesomie.GetQext(), mesomie.GetQsca()


def eval_mie_spectrum(x, m):
    Qext = []
    Qsca = []
    for i in range(len(x)):
        ext, sca =  eval_mie_single_spectrum_point(x[i], m[i])
        Qext.append(ext)
        Qsca.append(sca) 
    return Qext, Qsca

def eval_mesomie_spectrum(x, m, R, d_perp, d_parl=None, n_m=1):
    if np.array(d_perp == None).any():
        d_perp = np.zeros(len(x))
    if np.array(d_parl == None).any():
        d_parl = np.zeros(len(x))
    Qext = []
    Qsca = []
    for i in range(len(x)):
        ext, sca =  eval_mesomie_single_spectrum_point(x[i], m[i], n_m, R,
                                                d_perp[i], d_parl[i])
        Qext.append(ext)
        Qsca.append(sca) 
    return Qext, Qsca


