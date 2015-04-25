#! /bin/sh
#
#    Copyright (C) 2009-2015 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2015 Konstantin Ladutenko <kostyfisik@gmail.com>
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
#    using it, cite the following reference:
#    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
#        a multilayered sphere," Computer Physics Communications,
#        vol. 180, Nov. 2009, pp. 2348-2354.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This test case calculates the electric field in the 
# XY plane, for a silver nanoshell embedded in water.

# Refractive index values correspond to the wavelength
# where maximum of the surface plasmon resonance (and,
# hence, of electric field) is expected.

PROGRAM='../../../fieldnlay'

$PROGRAM -l 2 0.38989409 1.16177963 0.00000000 0.46787291 0.42850284 5.47718289 -p -1.55957635 1.55957635 501 -1.55957635 1.55957635 501 0.00000000 0.00000000 1
