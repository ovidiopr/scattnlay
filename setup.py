#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2019 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2019 Konstantin Ladutenko <kostyfisik@gmail.com>
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

__version__ = "2.4"
__title__ = "Calculation of the scattering of EM radiation by a multilayered sphere"
__mod__ = "python-scattnlay"
__author__ = "Ovidio Peña Rodríguez"
__email__ = "ovidio@bytesfall.com"
__url__ = "https://github.com/ovidiopr/scattnlay"
__download_url__ = (
    "https://github.com/ovidiopr/scattnlay/archive/v" + __version__ + ".0.tar.gz"
)

import os
import sys
import tempfile
from setuptools import setup, Extension
import numpy as np
import pybind11 as pb

# Add brew paths
extra_include_dirs = []
extra_library_dirs = []
if sys.platform == 'darwin':
    # Check for homebrew
    for prefix in ['/opt/homebrew', '/usr/local']:
        if os.path.isdir(prefix):
            inc = os.path.join(prefix, 'include')
            lib = os.path.join(prefix, 'lib')
            if os.path.isdir(inc): extra_include_dirs.append(inc)
            if os.path.isdir(lib): extra_library_dirs.append(lib)
            
            # Specific highway path if needed
            hwy_prefix = os.path.join(prefix, 'opt', 'highway')
            if os.path.isdir(hwy_prefix):
                extra_include_dirs.append(os.path.join(hwy_prefix, 'include'))
                extra_library_dirs.append(os.path.join(hwy_prefix, 'lib'))

# Helper to check for headers/libraries
def check_compilation(code, extra_args=None, libraries=None):
    try:
        import distutils.ccompiler
        import distutils.sysconfig
        from distutils.errors import CompileError, LinkError
    except ImportError:
        try:
            from setuptools._distutils import ccompiler, sysconfig
            from setuptools._distutils.errors import CompileError, LinkError
            distutils = type('distutils', (), {})
            distutils.ccompiler = ccompiler
            distutils.sysconfig = sysconfig
        except ImportError:
            print("Warning: Could not import distutils to check dependencies. Assuming missing.")
            return False

    compiler = distutils.ccompiler.new_compiler()
    distutils.sysconfig.customize_compiler(compiler)
    
    # Add include dirs (including current directory for relative includes)
    include_dirs = [np.get_include(), pb.get_include(), '.'] + extra_include_dirs
    
    with tempfile.TemporaryDirectory() as tmpdir:
        fname = os.path.join(tmpdir, 'test.cpp')
        with open(fname, 'w') as f:
            f.write(code)
        
        try:
            # Compile
            # Suppress output
            obj_files = compiler.compile([fname], output_dir=tmpdir, include_dirs=include_dirs, extra_postargs=extra_args)
            # Link (if libraries are needed)
            if libraries:
                compiler.link_executable(obj_files, os.path.join(tmpdir, 'test'), libraries=libraries, library_dirs=extra_library_dirs)
            return True
        except (CompileError, LinkError, Exception):
            return False

# Define extensions
ext_dp = Extension(
    "scattnlay_dp",
    ["src/pb11-wrapper.cc"],
    language="c++",
    include_dirs=[np.get_include(), pb.get_include(), '.'] + extra_include_dirs,
    library_dirs=extra_library_dirs,
    extra_compile_args=["-std=c++11"],
)

ext_mp = Extension(
    "scattnlay_mp",
    ["src/pb11-wrapper.cc"],
    language="c++",
    include_dirs=[np.get_include(), pb.get_include(), '.'] + extra_include_dirs,
    library_dirs=extra_library_dirs,
    extra_compile_args=["-std=c++14", "-DMULTI_PRECISION=100"],
)

ext_simd = Extension(
    "scattnlay_simd",
    ["src/pb11-wrapper.cc"],
    language="c++",
    include_dirs=[np.get_include(), pb.get_include(), '.'] + extra_include_dirs,
    library_dirs=extra_library_dirs,
    extra_compile_args=["-std=c++17", "-DWITH_HWY", "-DHWY_DISABLED_TARGETS=0x4000000"],
    libraries=["hwy"],
)

# Determine which extensions to build
extensions = [ext_dp]

# Check for Boost (MP)
print("Checking for Boost Multiprecision...")
boost_code = "#include <boost/multiprecision/cpp_bin_float.hpp>\nint main() { return 0; }"
if check_compilation(boost_code, extra_args=["-std=c++14"]):
    print("Boost found. Enabling scattnlay_mp.")
    extensions.append(ext_mp)
else:
    print("Boost not found. Skipping scattnlay_mp.")

# Check for Highway (SIMD)
print("Checking for Google Highway...")
hwy_code = "#include <hwy/highway.h>\nint main() { return 0; }"
if check_compilation(hwy_code, extra_args=["-std=c++17"], libraries=["hwy"]):
    print("Highway found. Enabling scattnlay_simd.")
    extensions.append(ext_simd)
else:
    print("Highway not found. Skipping scattnlay_simd.")

setup(
    name=__mod__,
    version=__version__,
    description=__title__,
    long_description="""The Python version of scattnlay, a computer implementation of the algorithm for the     calculation of electromagnetic radiation scattering by a multilayered sphere developed by Yang. It has been     shown that the program is effective, resulting in very accurate values of scattering efficiencies for a wide     range of size parameters, which is a considerable improvement over previous implementations of similar algorithms.     For details see: O. Pena, U. Pal, Comput. Phys. Commun. 180 (2009) 2348-2354.""",
    author=__author__,
    author_email=__email__,
    maintainer=__author__,
    maintainer_email=__email__,
    keywords=[
        "Mie scattering",
        "Multilayered sphere",
        "Efficiency factors",
        "Cross-sections",
    ],
    url=__url__,
    download_url=__download_url__,
    license="GPL",
    platforms="any",
    packages=["scattnlay"],
    ext_modules=extensions,
    install_requires=["numpy"],
)
