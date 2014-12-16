//**********************************************************************************//
//    Copyright (C) 2009-2015  Ovidio Pena <ovidio@bytesfall.com>                   //
//                                                                                  //
//    This file is part of scattnlay                                                //
//                                                                                  //
//    This program is free software: you can redistribute it and/or modify          //
//    it under the terms of the GNU General Public License as published by          //
//    the Free Software Foundation, either version 3 of the License, or             //
//    (at your option) any later version.                                           //
//                                                                                  //
//    This program is distributed in the hope that it will be useful,               //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 //
//    GNU General Public License for more details.                                  //
//                                                                                  //
//    The only additional remark is that we expect that all publications            //
//    describing work using this software, or all commercial products               //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//                                                                                  //
//    You should have received a copy of the GNU General Public License             //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//**********************************************************************************//

#define VERSION "0.3.1"
#include <complex>
#include <vector>

int ScattCoeffs(int L, int pl, std::vector<double> x, std::vector<std::complex<double> > m, int n_max,
		        std::vector<std::complex<double> > &an, std::vector<std::complex<double> > &bn);

int nMie(int L, std::vector<double> x, std::vector<std::complex<double> > m,
         int nTheta, std::vector<double> Theta,
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
         std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2);

int nMie(int L, int pl, std::vector<double> x, std::vector<std::complex<double> > m,
         int nTheta, std::vector<double> Theta,
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
         std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2);

int nMie(int L, std::vector<double> x, std::vector<std::complex<double> > m,
         int nTheta, std::vector<double> Theta, int n_max,
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
         std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2);

int nMie(int L, int pl, std::vector<double> x, std::vector<std::complex<double> > m,
         int nTheta, std::vector<double> Theta, int n_max,
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
		 std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2);

int nField(int L, int pl, std::vector<double> x, std::vector<std::complex<double> > m, int nmax,
           int ncoord, std::vector<double> Xp, std::vector<double> Yp, std::vector<double> Zp,
		   std::vector<std::vector<std::complex<double> > >& E, std::vector<std::vector<std::complex<double> > >& H);


