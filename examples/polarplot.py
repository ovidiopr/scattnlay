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

from scattnlay import scattnlay
import numpy as np
import cmath
import matplotlib.pyplot as plt
import matplotlib as mp
import math as m
import sys,os,shutil,glob,subprocess

def polarplot(x, m, fname, plot_type="unpolarized", mode_n = -1, mode_type = -1):
    theta = np.linspace(0,2*np.pi, 360*4)
    theta = theta[:-1]
    terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(
            np.array([x]), np.array([m]), theta =  theta, mode_n = mode_n, mode_type = mode_type)
    S = 1.0/2.0*( np.abs(S1[0])**2 + np.abs(S2[0])**2 )
    print(S)
    font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : '20'}
    mp.rc('font', **font)
    mp.rcParams['axes.linewidth'] = 2
    #mp.rcParams['grid.linewidth'] = 1.5
    mp.rcParams['legend.fontsize'] = 20
    mp.rcParams['axes.labelweight'] = 'bold'
    #mp.rcParams['lines.linewidth'] = 1.5

    ax = plt.subplot(111, polar=True)
    
    # low_lim = np.min(S)
    # hi_lim = np.max(S)
    # ax.set_ylim(low_lim, hi_lim)
    # yticks = np.arange(low_lim,hi_lim, 5)

    # ax.set_yticks(yticks)
    # ax.set_yticklabels(yticks,fontsize='18')
    for lb in ax.get_yticklabels():
        lb.set_fontsize(15)
    
    # # ax.plot(theta, r, color='r', linewidth=3)
    # # ax.set_rmax(2.0)
    ax.grid(linestyle="--")

    # ax.set_ylabel("dBsm",rotation='horizontal');
    #ax.yaxis.set_label_coords(1.05, .7)
    # # ax.set_title("dBsm", loc=[0.78,0.97])
    # # tick locations
    # thetaticks = np.arange(0,360,45)
    # # set ticklabels location at 1.3 times the axes' radius
    # ax.set_thetagrids(thetaticks)# The frac parameter was deprecated in version
    #                              # 2.1. Use tick padding via Axes.tick_params
    #                              # instead. , frac=1.15)
    # ax.xaxis.label.set_visible(False)
    # ax.tick_params(pad=15)
    # plt.setp(ax.get_xticklabels(), visible=False)
    # plt.setp(ax.get_yticklabels(), visible=False)
    #plt.setp(ax.get_xticklines()[-2:], visible=False)
    #ax.grid()
    # Init line data structure, will be rewritten before real plotting.
    line_tscs, = ax.plot(theta, np.log(S), 'r', linewidth=3, label=r"$S_{11}$")
    # handles, labels = ax.get_legend_handles_labels()
    # #lg = ax.legend(handles, labels, loc=[0.74,0.97])
    # #lg = ax.legend(handles, labels, loc=[-0.24,-0.12])
    # lg = ax.legend(handles, labels, loc=[-0.36,0.9])
    # lg.draw_frame(False)

    # #two curves per figure
    # line_bar.set_ydata(data_dBsm[i])
    plt.savefig(fname) #"RCS-"+description[i]+"-"+description[i+num_of_plots]+".svg", format="svg")
    plt.clf()
    plt.close()
