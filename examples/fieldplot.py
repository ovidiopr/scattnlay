#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2015 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2015  Konstantin Ladutenko <kostyfisik@gmail.com>
#
#    This file is part of python-scattnlay
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    The only additional remark is that we expect that all publications
#    describing work using this software, or all commercial products
#    using it, cite the following reference:
#    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Several functions to plot field and streamlines (power flow lines).

from scattnlay import fieldnlay, scattnlay
import numpy as np


###############################################################################
def GetCoords(crossplane, npts, factor, x):
    """
    crossplane: XZ, YZ, XY, or XYZ (half is XZ, half is YZ)
    npts: number of point in each direction
    factor: ratio of plotting size to outer size of the particle
    x: size parameters for particle layers
    """
    scan = np.linspace(-factor*x[-1], factor*x[-1], npts)
    zero = np.zeros(npts*npts, dtype = np.float64)

    if crossplane=='XZ':
        coordX, coordZ = np.meshgrid(scan, scan)
        coordX.resize(npts*npts)
        coordZ.resize(npts*npts)
        coordY = zero
    elif crossplane == 'YZ':
        coordY, coordZ = np.meshgrid(scan, scan)
        coordY.resize(npts*npts)
        coordZ.resize(npts*npts)
        coordX = zero
    elif crossplane == 'XY':
        coordX, coordY = np.meshgrid(scan, scan)
        coordX.resize(npts*npts)
        coordY.resize(npts*npts)
        coordZ = zero
    elif crossplane=='XYZ': # Upper half: XZ; Lower half: YZ
        coordX, coordZ = np.meshgrid(scan, scan)
        coordY, coordZ = np.meshgrid(scan, scan)
        coordX[:, scan<0] = 0
        coordY[:, scan>=0] = 0
        coordX.resize(npts*npts)
        coordY.resize(npts*npts)
        coordZ.resize(npts*npts)
    else:
        print("Unknown crossplane")
        import sys
        sys.exit()

    return coordX, coordY, coordZ, scan


###############################################################################
def GetField(crossplane, npts, factor, x, m, pl):
    """
    crossplane: XZ, YZ, XY, or XYZ (half is XZ, half is YZ)
    npts: number of point in each direction
    factor: ratio of plotting size to outer size of the particle
    x: size parameters for particle layers
    m: relative index values for particle layers
    """
    coordX, coordY, coordZ, scan = GetCoords(crossplane, npts, factor, x)

    terms, E, H = fieldnlay(x, m, coordX, coordY, coordZ, pl=pl)
    if len(E.shape) > 2:
        E = E[0, :, :]
        H = H[0, :, :]

    S = np.cross(E, np.conjugate(H)).real
    print(S)

    if crossplane=='XZ':
        Sx = np.resize(S[:, 2], (npts, npts)).T
        Sy = np.resize(S[:, 0], (npts, npts)).T
    elif crossplane == 'YZ':
        Sx = np.resize(S[:, 2], (npts, npts)).T
        Sy = np.resize(S[:, 1], (npts, npts)).T
    elif crossplane == 'XY':
        Sx = np.resize(S[:, 1], (npts, npts)).T
        Sy = np.resize(S[:, 0], (npts, npts)).T
    elif crossplane=='XYZ': # Upper half: XZ; Lower half: YZ
        Sx = np.resize(S[:, 2], (npts, npts)).T
        Sy = np.resize(S[:, 0], (npts, npts)).T
        Sy[scan<0] = np.resize(S[:, 1], (npts, npts)).T[scan<0]
    else:
        print("Unknown crossplane")
        import sys
        sys.exit()

    return E, H, S, scan, Sx, Sy
###############################################################################


def fieldplot(fig, ax, x, m, WL, comment='', WL_units=' ', crossplane='XZ',
              field_to_plot='Pabs', npts=101, factor=2.1,
              flow_total=11, density=20, minlength=0.1, maxlength=4.0,
              arrowstyle='-|>', arrowsize=1.0,
              pl=-1, draw_shell=False, outline_width=1, subplot_label=' '):

    E, H, S, scan, Sx, Sy = GetField(crossplane, npts, factor, x, m, pl)
    Er = np.absolute(E)
    Hr = np.absolute(H)
    try:
        from matplotlib import cm
        from matplotlib.colors import LogNorm

        if field_to_plot == 'Pabs':
            label = r'$\operatorname{Re}(E \times H^*)$'
            data = np.resize(np.linalg.norm(np.cross(E, np.conjugate(H)), axis=1).real, (npts, npts)).T
        elif field_to_plot == 'Eabs':
            label = r'$|E|$'
            Eabs = np.sqrt(Er[:, 0]**2 + Er[:, 1]**2 + Er[:, 2]**2)
            data = np.resize(Eabs, (npts, npts)).T
        elif field_to_plot == 'Habs':
            label = r'$|H|$'
            Habs = np.sqrt(Hr[:, 0]**2 + Hr[:, 1]**2 + Hr[:, 2]**2)
            Habs = 376.730313667*Habs # scale to free space impedance
            data = np.resize(Habs, (npts, npts)).T
        elif field_to_plot == 'angleEx':
            label = r'$arg(E_x)$'
            Eangle = np.angle(E[:, 0])/np.pi*180
            data = np.resize(Eangle, (npts, npts)).T
        elif field_to_plot == 'angleHy':
            label = r'$arg(H_y)$'
            Hangle = np.angle(H[:, 1])/np.pi*180
            data = np.resize(Hangle, (npts, npts)).T

        # Rescale to better show the axes
        scale = scan*WL/2.0/np.pi

        # Define scale ticks
        min_tick = np.amin(data[~np.isnan(data)])
        max_tick = np.amax(data[~np.isnan(data)])

        scale_ticks = np.linspace(min_tick, max_tick, 5)

        ax.set_title(label)
        # build a rectangle in axes coords
        ax.annotate(subplot_label, xy=(0.0, 1.1), xycoords='axes fraction',  # fontsize=10,
                    horizontalalignment='left', verticalalignment='top')

        # Interpolation can be 'nearest', 'bilinear' or 'bicubic'
        cax = ax.imshow(data, interpolation='nearest', cmap=cm.jet,
                        origin='lower', vmin=min_tick, vmax=max_tick,
                        extent=(min(scale), max(scale), min(scale), max(scale))
                        # ,norm = LogNorm()
                        )
        ax.axis("image")

        # Add colorbar
        cbar = fig.colorbar(cax, ticks=[a for a in scale_ticks], ax=ax)
        # vertically oriented colorbar
        if 'angle' in field_to_plot:
            cbar.ax.set_yticklabels(['%3.0f' % (a) for a in scale_ticks])
        else:
            cbar.ax.set_yticklabels(['%g' % (a) for a in scale_ticks])

        if crossplane == 'XZ':
            ax.set_xlabel('Z (%s)' % (WL_units))
            ax.set_ylabel('X (%s)' % (WL_units))
        elif crossplane == 'YZ':
            ax.set_xlabel('Z (%s)' % (WL_units))
            ax.set_ylabel('Y (%s)' % (WL_units))
        elif crossplane=='XYZ':
            ax.set_xlabel('Z (%s)' % (WL_units))
            ax.set_ylabel('Y(<0):X(>0) (%s)' % (WL_units))

            # draw a line to separate both planes
            ax.axhline(linewidth=outline_width, color='black')
        elif crossplane == 'XY':
            ax.set_xlabel('X (%s)' % (WL_units))
            ax.set_ylabel('Y (%s)' % (WL_units))

        if draw_shell:
            # Draw the nanoshell
            from matplotlib import patches
            from matplotlib.path import Path
            for xx in x:
                r = xx*WL/2.0/np.pi
                s1 = patches.Arc((0, 0), 2.0*r, 2.0*r,  angle=0.0, zorder=1.8,
                                 theta1=0.0, theta2=360.0, linewidth=outline_width, color='black')
                ax.add_patch(s1)

        # Draw flow lines
        if (not crossplane == 'XY') and flow_total > 0:
            margin = 0.98
            points = np.vstack((margin*scale.min()*np.ones(flow_total),
                                np.linspace(margin*scale.min(),
                                            margin*scale.max(), flow_total))).transpose()

            # Plot the streamlines with an appropriate colormap and arrow style
            ax.streamplot(scale, scale, Sx, Sy,
                          start_points=points, integration_direction='both',
                          density=density, minlength=minlength, maxlength=maxlength,
                          linewidth=outline_width, color='white',
                          arrowstyle=arrowstyle, arrowsize=arrowsize)
    finally:
        terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(x, m)
        print("Qsca = " + str(Qsca))
    #
