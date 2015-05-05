#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2015 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
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

import scattnlay
from scattnlay import fieldnlay
from scattnlay import scattnlay
import numpy as np
import cmath


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle
###############################################################################
def GetFlow3D(x0, y0, z0, max_length, max_angle, x, m):
    # Initial position
    flow_x = [x0]
    flow_y = [y0]
    flow_z = [z0]
    max_step = x[-1]/3
    min_step = x[0]/2000
#    max_step = min_step
    step = min_step*2.0
    if max_step < min_step:
        max_step = min_step
    coord = np.vstack(([flow_x[-1]], [flow_y[-1]], [flow_z[-1]])).transpose()
    terms, E, H = fieldnlay(np.array([x]), np.array([m]), coord)
    Ec, Hc = E[0, 0, :], H[0, 0, :]
    S = np.cross(Ec, Hc.conjugate()).real
    Snorm_prev = S/np.linalg.norm(S)
    length = 0
    dpos = step
    while length < max_length:
        if step<max_step:
                step = step*2.0
        r = np.sqrt(flow_x[-1]**2 + flow_y[-1]**2 + flow_z[-1]**2)
        while step > min_step:
            #Evaluate displacement from previous poynting vector
            dpos = step
            dx = dpos*Snorm_prev[0];
            dy = dpos*Snorm_prev[1];
            dz = dpos*Snorm_prev[2];
            #Test the next position not to turn more than max_angle
            coord = np.vstack(([flow_x[-1]+dx], [flow_y[-1]+dy], [flow_z[-1]+dz])).transpose()
            terms, E, H = fieldnlay(np.array([x]), np.array([m]), coord)
            Ec, Hc = E[0, 0, :], H[0, 0, :]
            Eth = max(np.absolute(Ec))/1e10
            Hth = max(np.absolute(Hc))/1e10
            for i in xrange(0,len(Ec)):
                if abs(Ec[i]) < Eth:
                    Ec[i] = 0+0j
                if abs(Hc[i]) < Hth:
                    Hc[i] = 0+0j
            S = np.cross(Ec, Hc.conjugate()).real
            Snorm = S/np.linalg.norm(S)
            # diff = Snorm-Snorm_prev
            # if np.linalg.norm(diff)<0.05:
            angle = angle_between(Snorm, Snorm_prev)
            if abs(angle) < max_angle:
                break
            step = step/2.0
        #3. Save result
        Snorm_prev = Snorm
        dx = dpos*Snorm_prev[0];
        dy = dpos*Snorm_prev[1];
        dz = dpos*Snorm_prev[2];
        length = length + step
        flow_x.append(flow_x[-1] + dx)
        flow_y.append(flow_y[-1] + dy)
        flow_z.append(flow_z[-1] + dz)

    return np.array(flow_x), np.array(flow_y), np.array(flow_z)


###############################################################################
def GetField(crossplane, npts, factor, x, m):
    """
    crossplane: XZ, YZ, XY
    npts: number of point in each direction
    factor: ratio of plotting size to outer size of the particle
    x: size parameters for particle layers
    m: relative index values for particle layers
    """
    scan = np.linspace(-factor*x[-1], factor*x[-1], npts)
    zero = np.zeros(npts*npts, dtype = np.float64)

    if crossplane=='XZ':
        coordX, coordZ = np.meshgrid(scan, scan)
        coordX.resize(npts*npts)
        coordZ.resize(npts*npts)
        coordY = zero
        coordPlot1 = coordX
        coordPlot2 = coordZ
    elif crossplane=='YZ':
        coordY, coordZ = np.meshgrid(scan, scan)
        coordY.resize(npts*npts)
        coordZ.resize(npts*npts)
        coordX = zero
        coordPlot1 = coordY
        coordPlot2 = coordZ
    elif crossplane=='XY':
        coordX, coordY = np.meshgrid(scan, scan)
        coordX.resize(npts*npts)
        coordY.resize(npts*npts)
        coordZ = zero
        coordPlot1 = coordY
        coordPlot2 = coordX
        
    coord = np.vstack((coordX, coordY, coordZ)).transpose()
    terms, E, H = fieldnlay(np.array([x]), np.array([m]), coord)
    Ec = E[0, :, :]
    Hc = H[0, :, :]
    P=[]
    P = np.array(map(lambda n: np.linalg.norm(np.cross(Ec[n], Hc[n])).real, range(0, len(E[0]))))

    # for n in range(0, len(E[0])):
    #     P.append(np.linalg.norm( np.cross(Ec[n], np.conjugate(Hc[n]) ).real/2 ))
    return Ec, Hc, P, coordPlot1, coordPlot2


