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
# E-k plane, for an spherical Si-Ag-Si nanoparticle. Core radius is 17.74 nm,
# inner layer 23.31nm, outer layer 22.95nm. Working wavelength is 800nm, we use
# silicon epsilon=13.64+i0.047, silver epsilon= -28.05+i1.525

import os, cmath
import numpy as np
from scattnlay import fieldnlay
from fieldplot import fieldplot


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("dirnames", nargs='*', default='.', help="read all data from DIR(S)")
    parser.add_argument("-f", "--filename", dest="fname", nargs='?', default=None,
                        help="name of 'n' file")
    parser.add_argument("-w", "--wavelength", dest="wl", default=3.75, type=float,
                        help="wavelength of electromagnetic wave")
    parser.add_argument("-r", "--radius", dest="rad", default=None, type=float,
                        help="radius of PEC sphere")
    parser.add_argument("-t", "--thickness", dest="tc", default=0.8, type=float,
                        help="thickness of cloaking layer")
    parser.add_argument("-n", "--npoints", dest="npts", default=101, type=int,
                        help="number of points for the grid")

    args = parser.parse_args()

    for dirname in args.dirnames:
        print "Calculating spectra for data file(s) in dir '%s'..." % (dirname)

        wl = args.wl # cm
        if (args.rad is None):
            Rs = 0.75*wl  # cm
        else:
            Rs = args.rad # cm
        tc = args.tc # cm

        if (args.fname is None):
            files = [x for x in os.listdir('%s/' % (dirname)) if x.endswith('.dat')]
            files.sort()
        else:
            files = [args.fname]

        npts = args.npts # cm

        if not os.path.exists('%s/flow-results/' % (dirname)):
            os.makedirs('%s/flow-results/' % (dirname))

        Rt = Rs + tc # cm

        print "Wl = %.2f, Rs = %.2f, tc = %.2f, Rt = %.2f" % (wl, Rs, tc, Rt)

        ms = 1.0 + 40.0j
        for i, fname in enumerate(files):
            print "Calculating spectra for file '%s'..." % (fname)

            basename = os.path.splitext(fname)[0]

            nvalues = np.loadtxt('%s/%s' % (dirname, fname))*1.0 + 1e-11j

            tl = tc/len(nvalues) # cm
            r = [Rs]

            for i in range(len(nvalues)):
                r += [r[i] + tl]

            x = np.ones((len(nvalues) + 1), dtype = np.float64)
            m = np.ones((len(nvalues) + 1), dtype = np.complex128)

            x = 2.0*np.pi*np.array(r, dtype = np.float64)/wl
            m = np.array([ms] + nvalues[:, 1].tolist(), dtype = np.complex128)
            print(x,m)

            factor = 2.91*x[0]/x[-1]
            print factor
            comment='PEC-'+basename
            WL_units=''
            #flow_total = 39
            # flow_total = 23 #SV False
            flow_total = 24
            #flow_total = 4
            #crossplane='XZ'
            crossplane='XYZ'
            #crossplane='YZ'
            #crossplane='XY'

            # Options to plot: Eabs, Habs, Pabs, angleEx, angleHy
            #field_to_plot='Pabs'
            #field_to_plot='Eabs'
            
            field_to_plot='angleEx'
            #field_to_plot='angleHy'
            print "x =", x
            print "m =", m

            import matplotlib.pyplot as plt
            plt.rcParams.update({'font.size': 16})
            fig, axs = plt.subplots(1,1)#, sharey=True, sharex=True)
            fig.tight_layout()
            fieldplot(fig, axs, x,m, wl, comment, WL_units, crossplane, field_to_plot, npts, factor, flow_total,
                      subplot_label=' ',is_flow_extend=False
                      , outline_width=1.5
                      , pl=0 #PEC layer starts the design
                      )
            # fieldplot(fig, axs, x[0],m[0], wl, comment, WL_units, crossplane, field_to_plot, npts, factor, flow_total,
            #           subplot_label=' ' ,is_flow_extend=False
            #           , outline_width=1.5
            #           , pl=0 #PEC layer starts the design
            #           )
            fig.subplots_adjust(hspace=0.3, wspace=-0.1)
            plt.savefig(comment+"-R"+str(int(round(x[-1]*wl/2.0/np.pi)))+"-"+crossplane+"-"
                        +field_to_plot+".pdf",pad_inches=0.02, bbox_inches='tight')
            plt.draw()
            plt.clf()
            plt.close()


        print "Done!!"

