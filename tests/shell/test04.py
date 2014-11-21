#! /usr/bin/python
# This python script calculates the differential cross section
# from the  scattering amplitudes values given by scattnlay.
# It expects the values to be in a file named 'luneburg.txt'
# and the output will be saved to a file named 'luneburg_fin.txt'

import scipy
import pylab

def diff_scattering(Theta):
  return scipy.cos(Theta*scipy.pi/180.0)

data = scipy.loadtxt(fname = 'luneburg.txt', delimiter = ', ', skiprows = 2)

S11 = data[ : , 1]*data[ : , 1] + data[ : , 2]*data[ : , 2] + data[ : , 3]*data[ : , 3] + data[ : , 4]*data[ : , 4]
S11 = S11/7200.0
S11 = scipy.vstack((data[ : , 0], S11, diff_scattering(data[ : , 0]))).transpose()

scipy.savetxt('luneburg_fin.txt', S11)

pylab.plot(S11[ : , 0], S11[ : , 1], S11[ : , 0], S11[ : , 2], color='k')

ax = pylab.gca()
ax.set_yscale('log')
ax.set_ylim(1e-4, 1e3)

pylab.draw()

pylab.show()


