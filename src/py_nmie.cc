//**********************************************************************************//
//    Copyright (C) 2009-2018  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2018  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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
//    using it, cite at least one of the following references:                      //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie              //
//        calculation of electromagnetic near-field for a multilayered              //
//        sphere," Computer Physics Communications, vol. 214, May 2017,             //
//        pp. 225-230.                                                              //
//                                                                                  //
//    You should have received a copy of the GNU General Public License             //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//**********************************************************************************//

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "nmie.hpp"
#include "py_nmie.h"

// Same as ScattCoeffs in 'nmie.h' but uses double arrays to return the results (useful for python).
// This is a workaround because I have not been able to return the results using 
// std::vector<std::complex<double> >
int ScattCoeffs(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m,
                const int nmax, double anr[], double ani[], double bnr[], double bni[]) {

  int i, result;
  std::vector<std::complex<double> > an, bn;
  an.resize(nmax);
  bn.resize(nmax);

  result = nmie::ScattCoeffs(L, pl, x, m, nmax, an, bn);

  for (i = 0; i < result; i++) {
    anr[i] = an[i].real();
    ani[i] = an[i].imag();
    bnr[i] = bn[i].real();
    bni[i] = bn[i].imag();
  }

  return result;
}

// Same as nMie in 'nmie.h' but uses double arrays to return the results (useful for python).
// This is a workaround because I have not been able to return the results using 
// std::vector<std::complex<double> >
int nMie(const int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m,
         const int nTheta, std::vector<double>& Theta, const int nmax,
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
		 double S1r[], double S1i[], double S2r[], double S2i[]) {

  int i, result;
  std::vector<std::complex<double> > S1, S2;
  S1.resize(nTheta);
  S2.resize(nTheta);

  result = nmie::nMie(L, pl, x, m, nTheta, Theta, nmax, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);

  for (i = 0; i < nTheta; i++) {
    S1r[i] = S1[i].real();
    S1i[i] = S1[i].imag();
    S2r[i] = S2[i].real();
    S2i[i] = S2[i].imag();
  }

  return result;
}

// Same as nField in 'nmie.h' but uses double arrays to return the results (useful for python).
// This is a workaround because I have not been able to return the results using 
// std::vector<std::complex<double> >
int nField(const int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const int nmax,
           const int nCoords, std::vector<double>& Xp, std::vector<double>& Yp, std::vector<double>& Zp,
           double Erx[], double Ery[], double Erz[], double Eix[], double Eiy[], double Eiz[],
           double Hrx[], double Hry[], double Hrz[], double Hix[], double Hiy[], double Hiz[]) {

  int i, result;
  std::vector<std::vector<std::complex<double> > > E, H;
  E.resize(nCoords);
  H.resize(nCoords);
  for (i = 0; i < nCoords; i++) {
    E[i].resize(3);
    H[i].resize(3);
  }

  result = nmie::nField(L, pl, x, m, nmax, nCoords, Xp, Yp, Zp, E, H);

  for (i = 0; i < nCoords; i++) {
    Erx[i] = E[i][0].real();
    Ery[i] = E[i][1].real();
    Erz[i] = E[i][2].real();
    Eix[i] = E[i][0].imag();
    Eiy[i] = E[i][1].imag();
    Eiz[i] = E[i][2].imag();
    Hrx[i] = H[i][0].real();
    Hry[i] = H[i][1].real();
    Hrz[i] = H[i][2].real();
    Hix[i] = H[i][0].imag();
    Hiy[i] = H[i][1].imag();
    Hiz[i] = H[i][2].imag();
  }

  return result;
}

