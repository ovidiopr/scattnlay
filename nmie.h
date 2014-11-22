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
#include <string>
#include <vector>

int nMie(int L, double x[], complex m[], int nTheta, double Theta[], double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[]);

int nMiePEC(int L, int pl, double x[], complex m[], int nTheta, double Theta[], double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[]);

int nMieMax(int L, double x[], complex m[], int nTheta, double Theta[], int n_max, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[]);

int nMieScatt(int L, int pl, double x[], complex m[], int nTheta, double Theta[], int n_max, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[]);

int nMieField(int L, int pl, double x[], complex m[], int n_max, int nCoords, double Xp[], double Yp[], double Zp[], complex E[], complex H[]);


int nMie_std(double x[], complex m[], double Theta[], complex S1[], complex S2[],
int L, std::vector<double> &x_std, std::vector<std::complex<double> > &m_std, int nTheta, std::vector<double> &Theta_std, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector< std::complex<double> > &S1_std, std::vector< std::complex<double> >  &S2_std);


int nMieScatt_std(double x[], complex m[], double Theta[], complex S1[], complex S2[],
		  int L, int pl,
		  std::vector<double> &x_std, std::vector<std::complex<double> > &m_std,
		  int nTheta, std::vector<double> &Theta_std,
		  int n_max, double *Qext, double *Qsca, double *Qabs, double *Qbk,
		  double *Qpr, double *g, double *Albedo,
		  std::vector< std::complex<double> > &S1_std,
		  std::vector< std::complex<double> >  &S2_std);





