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

#include <stdlib.h>
#include <stdio.h>
#include "ucomplex.h"
#include "nmie.h"
#include "py_nmie.h"

//Same as nMieScatt but only uses doubles (useful for python)
int nfMieScatt(int L, int pl, double x[], double mr[], double mi[], int nTheta, double Theta[], int n_max, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, double S1r[], double S1i[], double S2r[], double S2i[]){
  int i, result;
  complex m[L], S1[nTheta], S2[nTheta];

  for (i = 0; i < L; i++) {
    m[i].r = mr[i];
    m[i].i = mi[i];
  }

  result = nMieScatt(L, pl, x, m, nTheta, Theta, n_max, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);

  for (i = 0; i < nTheta; i++) {
    S1r[i] = S1[i].r;
    S1i[i] = S1[i].i;
    S2r[i] = S2[i].r;
    S2i[i] = S2[i].i;
  }

  return result;
}

int nfMieField(int L, int pl, double x[], double mr[], double mi[], int n_max, int nCoords, double Xp[], double Yp[], double Zp[], double Er[], double Ei[], double Hr[], double Hi[]){
  int i, result;
  complex m[L], E[nCoords], H[nCoords];

  for (i = 0; i < L; i++) {
    m[i].r = mr[i];
    m[i].i = mi[i];
  }

  result = nMieField(L, pl, x, m, n_max, nCoords, Xp, Yp, Zp, E, H);

  for (i = 0; i < nCoords; i++) {
    Er[i] = E[i].r;
    Ei[i] = E[i].i;
    Hr[i] = H[i].r;
    Hi[i] = H[i].i;
  }

  return result;
}

