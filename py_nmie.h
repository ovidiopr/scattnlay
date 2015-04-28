//**********************************************************************************//
//    Copyright (C) 2009-2015  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2015  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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

#include <complex>
#include <vector>

using namespace std;

int ScattCoeffs(const unsigned int L, const int pl, vector<double>& x, vector<complex<double> >& m,
                const int nmax, double anr[], double ani[], double bnr[], double bni[]);

int nMie(const int L, const int pl, vector<double>& x, vector<complex<double> >& m,
         const int nTheta, vector<double>& Theta, const int nmax,
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
		 double S1r[], double S1i[], double S2r[], double S2i[]);

int nField(const int L, const int pl, vector<double>& x, vector<complex<double> >& m, const int nmax,
           const int nCoords, vector<double>& Xp, vector<double>& Yp, vector<double>& Zp,
           double Erx[], double Ery[], double Erz[], double Eix[], double Eiy[], double Eiz[],
           double Hrx[], double Hry[], double Hrz[], double Hix[], double Hiy[], double Hiz[]);

