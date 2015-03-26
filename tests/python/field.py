#!/usr/bin/env python

# This test case calculates the differential scattering
# cross section from a Luneburg lens, as described in:
# B. R. Johnson, Applied Optics 35 (1996) 3286-3296.

# The Luneburg lens is a sphere of radius a, with a
# radially-varying index of refraction, given by:
# m(r) = [2 - (r/a)**1]**(1/2)

# For the calculations, the Luneburg lens was approximated
# as a multilayered sphere with 500 equally spaced layers.
# The refractive index of each layer is defined to be equal to
# m(r) at the midpoint of the layer: ml = [2 - (xm/xL)**1]**(1/2),
# with xm = (xl-1 + xl)/2, for l = 1,2,...,L. The size
# parameter in the lth layer is xl = l*xL/500. According to
# geometrical optics theory, the differential cross section
# can be expressed as:
# d(Csca)/d(a**2*Omega) = cos(Theta)

# The differential cross section from wave optics is:
# d(Csca)/d(a**2*Omega) = S11(Theta)/x**2

from scattnlay import fieldnlay
import numpy as np

x = np.ones((1, 1), dtype = np.float64)
x[0, 0] = 0.1

m = np.ones((1, 1), dtype = np.complex128)
m[0, 0] = (0.05 + 2.070j)/1.46

npts = 1001

scan = np.linspace(-3.0*x[0, 0], 3.0*x[0, 0], npts)

coordX, coordY = np.meshgrid(scan, scan)
coordX.resize(npts*npts)
coordY.resize(npts*npts)
coordZ = np.zeros(npts*npts, dtype = np.float64)

coord = np.vstack((coordX, coordY, coordZ)).transpose()

terms, E, H = fieldnlay(x, m, coord)

Er = np.absolute(E)

# |E|/|Eo|
Eh = np.sqrt(Er[0, :, 0]**2 + Er[0, :, 1]**2 + Er[0, :, 2]**2)

result = np.vstack((coordX, coordY, coordZ, Eh)).transpose()

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import LogNorm

    min_tick = 0.1
    max_tick = 1.0

    edata = np.resize(Eh, (npts, npts))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Rescale to better show the axes
    scale_x = np.linspace(min(coordX), max(coordX), npts)
    scale_y = np.linspace(min(coordY), max(coordY), npts)

    # Define scale ticks
    min_tick = max(0.5, min(min_tick, np.amin(edata)))
    max_tick = max(max_tick, np.amax(edata))
    scale_ticks = np.power(10.0, np.linspace(np.log10(min_tick), np.log10(max_tick), 6))

    # Interpolation can be 'nearest', 'bilinear' or 'bicubic'
    cax = ax.imshow(edata, interpolation = 'bicubic', cmap = cm.afmhot,
                    origin = 'lower', vmin = min_tick, vmax = max_tick,
                    extent = (min(scale_x), max(scale_x), min(scale_y), max(scale_y)),
                    norm = LogNorm())

    # Add colorbar
    cbar = fig.colorbar(cax, ticks = [a for a in scale_ticks])
    cbar.ax.set_yticklabels(['%3.1e' % (a) for a in scale_ticks]) # vertically oriented colorbar
    pos = list(cbar.ax.get_position().bounds)
    fig.text(pos[0] - 0.02, 0.925, '|E|/|E$_0$|', fontsize = 14)

    plt.xlabel('X')
    plt.ylabel('Y')

    plt.draw()

    plt.show()

    plt.clf()
    plt.close()
finally:
    np.savetxt("field.txt", result, fmt = "%.5f")
    print result


