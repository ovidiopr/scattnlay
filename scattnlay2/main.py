from scattnlay2.pynmie import scattcoeffs_, scattnlay_,  fieldnlay_
import numpy as np
def scattcoeffs(x, m, nmax=-1, pl=-1, nstore=1000):
    if len(x.shape) != 2:
        raise ValueError('The size parameter (x) should be 2-D NumPy array.')
    if len(m.shape) != 2:
        raise ValueError('The relative refractive index (m) should be 2-D NumPy array.')
    if nmax != -1: nstore = nmax
    total_evaluations = len(x)
    terms = np.zeros((total_evaluations), dtype=int)
    an = np.zeros((total_evaluations, nstore), dtype =complex)
    bn = np.zeros((total_evaluations, nstore), dtype =complex)
    for i in range(total_evaluations):
        terms[i], a, b = scattcoeffs_(x[i], m[i], nmax=nmax, pl=pl)
        if terms[i] > nstore:
            raise ValueError('Not enougth storage was allocated to keep all terms.\n'
                             'Increase value of "nstore" argument in "scattcoeffs" function!')
        an[i,:terms[i]] = a
        bn[i,:terms[i]] = b
    return terms, an, bn
        
#scattcoeffs()
