#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2018 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2018 Konstantin Ladutenko <kostyfisik@gmail.com>
#
#    This file is part of scattnlay
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
#    using it, cite at least one of the following references:
#    [1] O. Peña and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Peña-Rodríguez, "Mie
#        calculation of electromagnetic near-field for a multilayered
#        sphere," Computer Physics Communications, vol. 214, May 2017,
#        pp. 225-230.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

__version__ = '2.2'
__title__ = 'Calculation of the scattering of EM radiation by a multilayered sphere'
__mod__ = 'python-scattnlay'
__author__ = 'Ovidio Peña Rodríguez'
__email__ = 'ovidio@bytesfall.com'
__url__ = 'https://github.com/ovidiopr/scattnlay'
__download_url__ = 'https://github.com/ovidiopr/scattnlay/archive/v2.2.0.tar.gz'

from distutils.core import setup
from distutils.extension import Extension
import numpy as np

setup(name = __mod__,
      version = __version__,
      description = __title__,
      long_description="""The Python version of scattnlay, a computer implementation of the algorithm for the calculation of electromagnetic \
radiation scattering by a multilayered sphere developed by Yang. It has been shown that the program is effective, \
resulting in very accurate values of scattering efficiencies for a wide range of size parameters, which is a \
considerable improvement over previous implementations of similar algorithms. For details see: \
O. Pena, U. Pal, Comput. Phys. Commun. 180 (2009) 2348-2354.""",
      author = __author__,
      author_email = __email__,
      maintainer = __author__,
      maintainer_email = __email__,
      keywords = ['Mie scattering', 'Multilayered sphere', 'Efficiency factors', 'Cross-sections'],
      url = __url__,
      download_url = __download_url__,
      license = 'GPL',
      platforms = 'any',
      ext_modules = [Extension("scattnlay",
                               ["src/nmie.cc", "src/py_nmie.cc", "src/scattnlay.cpp"],
                               language = "c++",
                               include_dirs = [np.get_include()], 
                               extra_compile_args=['-std=c++11']),
                     Extension("scattnlay_mp",
                               ["src/nmie.cc", "src/py_nmie.cc", "src/scattnlay_mp.cpp"],
                               language = "c++",
                               include_dirs = [np.get_include()], 
                               extra_compile_args=['-std=c++11', '-DMULTI_PRECISION=100'])]
)

