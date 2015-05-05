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

# This test case calculates the electric field in the 
# E-k plane, for an spherical Si-Ag-Si nanoparticle.

import scattnlay
from scattnlay import fieldnlay
from scattnlay import scattnlay
import numpy as np
import cmath
from fieldplot import GetFlow3D
from fieldplot import GetField

###############################################################################
def SetXM(design):
    """ design value:
    1: AgSi - a1
    2: SiAgSi - a1, b1
    3: SiAgSi - a1, b2
    """
    epsilon_Si = 18.4631066585 + 0.6259727805j
    epsilon_Ag = -8.5014154589 + 0.7585845411j
    index_Si = np.sqrt(epsilon_Si)
    index_Ag = np.sqrt(epsilon_Ag)
    isSiAgSi=True
    isBulk = False
    if design==1:
        #36	5.62055	0	31.93	4.06	49	5.62055	500
        isSiAgSi=False
        WL=500 #nm
        core_width = 0.0 #nm Si
        inner_width = 31.93 #nm Ag
        outer_width = 4.06 #nm  Si
    elif design==2:
        #62.5	4.48866	29.44	10.33	22.73	0	4.48866	500
        WL=500 #nm
        core_width = 29.44 #nm Si
        inner_width = 10.33 #nm Ag
        outer_width = 22.73 #nm  Si
    elif design == 3:
        #81.4	3.14156	5.27	8.22	67.91	0	3.14156	500
        WL=500 #nm
        core_width = 5.27 #nm Si
        inner_width = 8.22 #nm Ag
        outer_width = 67.91 #nm  Si
    elif design==4:
        WL=800 #nm
        epsilon_Si = 13.64 + 0.047j
        epsilon_Ag = -28.05 + 1.525j
        core_width = 17.74 #nm Si
        inner_width = 23.31 #nm Ag
        outer_width = 22.95 #nm  Si
    elif design==5:
        WL=354 #nm
        core_r = WL/20.0
        epsilon_Ag = -2.0 + 0.28j   #original
        index_Ag = np.sqrt(epsilon_Ag)
        x = np.ones((1), dtype = np.float64)
        x[0] = 2.0*np.pi*core_r/WL
        m = np.ones((1), dtype = np.complex128)
        m[0] = index_Ag
        # x = np.ones((2), dtype = np.float64)
        # x[0] = 2.0*np.pi*core_r/WL/4.0*3.0
        # x[1] = 2.0*np.pi*core_r/WL
        # m = np.ones((2), dtype = np.complex128)
        # m[0] = index_Ag
        # m[1] = index_Ag
        return x, m, WL


    core_r = core_width
    inner_r = core_r+inner_width
    outer_r = inner_r+outer_width

    nm = 1.0
    if isSiAgSi:
        x = np.ones((3), dtype = np.float64)
        x[0] = 2.0*np.pi*core_r/WL
        x[1] = 2.0*np.pi*inner_r/WL
        x[2] = 2.0*np.pi*outer_r/WL
        m = np.ones((3), dtype = np.complex128)
        m[0] = index_Si/nm
        m[1] = index_Ag/nm
    #    m[0, 1] = index_Si/nm
        m[2] = index_Si/nm
    else:
        # bilayer
        x = np.ones((2), dtype = np.float64)
        x[0] = 2.0*np.pi*inner_r/WL
        x[1] = 2.0*np.pi*outer_r/WL
        m = np.ones((2), dtype = np.complex128)
        m[0] = index_Ag/nm
        m[1] = index_Si/nm
    return x, m, WL


###############################################################################
#design = 1 #AgSi
#design = 2
#design = 3
#design = 4   # WL=800
design = 5   # Bulk Ag
x, m, WL = SetXM(design)


WL_units='nm'
comment='P-SiAgSi-flow'
comment='bulk-P-Ag-flow'
print "x =", x
print "m =", m
npts = 101
factor=2.2
flow_total = 3
#flow_total = 0
crossplane='XZ'
#crossplane='YZ'
#crossplane='XY'



Ec, Hc, P, coordX, coordZ = GetField(crossplane, npts, factor, x, m)

Er = np.absolute(Ec)
Hr = np.absolute(Hc)
# |E|/|Eo|
Eabs = np.sqrt(Er[ :, 0]**2 + Er[ :, 1]**2 + Er[ :, 2]**2)
Eangle = np.angle(Ec[ :, 0])/np.pi*180
Habs= np.sqrt(Hr[ :, 0]**2 + Hr[ :, 1]**2 + Hr[ :, 2]**2)
Hangle = np.angle(Hc[ :, 1])/np.pi*180

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import LogNorm

    Eabs_data = np.resize(P, (npts, npts)).T
    #Eabs_data = np.resize(Pabs, (npts, npts)).T
    # Eangle_data = np.resize(Eangle, (npts, npts)).T
    # Habs_data = np.resize(Habs, (npts, npts)).T
    # Hangle_data = np.resize(Hangle, (npts, npts)).T

    fig, ax = plt.subplots(1,1)
    # Rescale to better show the axes
    scale_x = np.linspace(min(coordX)*WL/2.0/np.pi, max(coordX)*WL/2.0/np.pi, npts)
    scale_z = np.linspace(min(coordZ)*WL/2.0/np.pi, max(coordZ)*WL/2.0/np.pi, npts)

    # Define scale ticks
    min_tick = np.amin(Eabs_data[~np.isnan(Eabs_data)])
    max_tick = np.amax(Eabs_data[~np.isnan(Eabs_data)])
    scale_ticks = np.linspace(min_tick, max_tick, 6)

    # Interpolation can be 'nearest', 'bilinear' or 'bicubic'
    ax.set_title('Pabs')
    cax = ax.imshow(Eabs_data, interpolation = 'nearest', cmap = cm.jet,
                    origin = 'lower'
                    , vmin = min_tick, vmax = max_tick
                    , extent = (min(scale_x), max(scale_x), min(scale_z), max(scale_z))
                    #,norm = LogNorm()
                    )
    ax.axis("image")

    # Add colorbar
    cbar = fig.colorbar(cax, ticks = [a for a in scale_ticks])
    cbar.ax.set_yticklabels(['%5.3g' % (a) for a in scale_ticks]) # vertically oriented colorbar
    pos = list(cbar.ax.get_position().bounds)
    #fig.text(pos[0] - 0.02, 0.925, '|E|/|E$_0$|', fontsize = 14)
    if crossplane=='XZ':
        plt.xlabel('Z, '+WL_units)
        plt.ylabel('X, '+WL_units)
    elif crossplane=='YZ':
        plt.xlabel('Z, '+WL_units)
        plt.ylabel('Y, '+WL_units)
    elif crossplane=='XY':
        plt.xlabel('Y, '+WL_units)
        plt.ylabel('X, '+WL_units)
    

    # # This part draws the nanoshell
    from matplotlib import patches
    from matplotlib.path import Path
    for xx in x:
        r= xx*WL/2.0/np.pi
        s1 = patches.Arc((0, 0), 2.0*r, 2.0*r,  angle=0.0, zorder=1.8,
                         theta1=0.0, theta2=360.0, linewidth=1, color='black')
        ax.add_patch(s1)
    # 
    # for flow in range(0,flow_total):
    #     flow_x, flow_z = GetFlow(scale_x, scale_z, Ec, Hc,
    #                              min(scale_x)+flow*(scale_x[-1]-scale_x[0])/(flow_total-1),
    #                              min(scale_z),
    #                              #0.0,
    #                              npts*16)
    #     verts = np.vstack((flow_z, flow_x)).transpose().tolist()
    #     #codes = [Path.CURVE4]*len(verts)
    #     codes = [Path.LINETO]*len(verts)
    #     codes[0] = Path.MOVETO
    #     path = Path(verts, codes)
    #     patch = patches.PathPatch(path, facecolor='none', lw=1, edgecolor='yellow')
    #     ax.add_patch(patch)
    if (crossplane=='XZ' or crossplane=='YZ') and flow_total>0:

        from matplotlib.path import Path
        scanSP = np.linspace(-factor*x[-1], factor*x[-1], npts)
        min_SP = -factor*x[-1]
        step_SP = 2.0*factor*x[-1]/(flow_total-1)
        x0, y0, z0 = 0, 0, 0
        max_length=factor*x[-1]*15
        #max_length=factor*x[-1]*5
        max_angle = np.pi/200
        #for flow in range(0,flow_total*2+1):
        for flow in range(0,flow_total):
            if crossplane=='XZ':
                #x0 = min_SP*2 + flow*step_SP
                x0 = min_SP + flow*step_SP
                z0 = min_SP
                #y0 = x[-1]/20 
            elif crossplane=='YZ':
                #y0 = min_SP*2 + flow*step_SP
                y0 = min_SP + flow*step_SP
                z0 = min_SP
                #x0 = x[-1]/20
            flow_xSP, flow_ySP, flow_zSP = GetFlow3D(x0, y0, z0, max_length, max_angle, x, m)
            if crossplane=='XZ':
                flow_z_plot = flow_zSP*WL/2.0/np.pi
                flow_f_plot = flow_xSP*WL/2.0/np.pi
            elif crossplane=='YZ':
                flow_z_plot = flow_zSP*WL/2.0/np.pi
                flow_f_plot = flow_ySP*WL/2.0/np.pi

            verts = np.vstack((flow_z_plot, flow_f_plot)).transpose().tolist()
            codes = [Path.LINETO]*len(verts)
            codes[0] = Path.MOVETO
            path = Path(verts, codes)
            #patch = patches.PathPatch(path, facecolor='none', lw=0.2, edgecolor='white',zorder = 2.7)
            patch = patches.PathPatch(path, facecolor='none', lw=1, edgecolor='white',zorder = 1.9)
            ax.add_patch(patch)
            ax.plot(flow_z_plot, flow_f_plot, 'x',ms=2, mew=0.1, linewidth=0.5, color='k', fillstyle='none')

    plt.savefig(comment+"-R"+str(int(round(x[-1]*WL/2.0/np.pi)))+"-"+crossplane+".svg")
    plt.draw()

#    plt.show()

    plt.clf()
    plt.close()
finally:
    terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(np.array([x]),
                                                                     np.array([m]))
    print("Qabs = "+str(Qabs));
#


