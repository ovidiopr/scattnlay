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

//**********************************************************************************//
// This library implements the algorithm for a multilayered sphere described by:    //
//    [1] W. Yang, "Improved recursive algorithm for light scattering by a          //
//        multilayered sphere,” Applied Optics,  vol. 42, Mar. 2003, pp. 1710-1720. //
//                                                                                  //
// You can find the description of all the used equations in:                       //
//    [2] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//                                                                                  //
// Hereinafter all equations numbers refer to [2]                                   //
//**********************************************************************************//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "nmie.h"

#define round(x) ((x) >= 0 ? (int)((x) + 0.5):(int)((x) - 0.5))

// Calculate Nstop - equation (17)
int Nstop(double xL) {
  int result;

  if (xL <= 8) {
    result = round(xL + 4*pow(xL, 1/3) + 1);
  } else if (xL <= 4200) {
    result = round(xL + 4.05*pow(xL, 1/3) + 2);
  } else {
    result = round(xL + 4*pow(xL, 1/3) + 2);
  }

  return result;
}

//**********************************************************************************//
int Nmax(int L, int fl, int pl,
         std::vector<double> x,
		 std::vector<std::complex<double> > m) {
  int i, result, ri, riM1;
  result = Nstop(x[L - 1]);
  for (i = fl; i < L; i++) {
    if (i > pl) {
      ri = round(std::abs(x[i]*m[i]));
    } else {
      ri = 0;
    }
    if (result < ri) {
      result = ri;
    }

    if ((i > fl) && ((i - 1) > pl)) {
      riM1 = round(std::abs(x[i - 1]* m[i]));
    } else {
      riM1 = 0;
    }
    if (result < riM1) {
      result = riM1;
    }
  }
  return result + 15;
}

//**********************************************************************************//
// This function calculates the spherical Bessel functions (jn and yn) for a given  //
// real value r. See pag. 87 B&H.                                                   //
//                                                                                  //
// Input parameters:                                                                //
//   r: Real argument to evaluate jn and yn                                         //
//   n_max: Maximum number of terms to calculate jn and yn                          //
//                                                                                  //
// Output parameters:                                                               //
//   jn, yn: Spherical Bessel functions (double)                                    //
//**********************************************************************************//
void sphericalBessel(double r, int n_max, std::vector<double> &j, std::vector<double> &y) {
  int n;

  if (n_max >= 1) {
    j[0] = sin(r)/r;
    y[0] = -cos(r)/r;
  }

  if (n_max >= 2) {
    j[1] = sin(r)/r/r - cos(r)/r;
    y[1] = -cos(r)/r/r - sin(r)/r;
  }

  for (n = 2; n < n_max; n++) {
    j[n] = double(n + n + 1)*j[n - 1]/r - j[n - 2];
    y[n] = double(n + n + 1)*y[n - 1]/r - h[n - 2];
  }
}

//**********************************************************************************//
// This function calculates the spherical Hankel functions (h1n and h2n) for a      //
// given real value r. See eqs. (4.13) and (4.14), pag. 87 B&H.                     //
//                                                                                  //
// Input parameters:                                                                //
//   r: Real argument to evaluate h1n and h2n                                       //
//   n_max: Maximum number of terms to calculate h1n and h2n                        //
//                                                                                  //
// Output parameters:                                                               //
//   h1n, h2n: Spherical Hankel functions (complex)                                 //
//**********************************************************************************//
void sphericalHankel(double r, int n_max, std::vector<std::complex<double> > &h1, std::vector<std::complex<double> > &h2) {
  int n;
  std::complex<double> j, y;
  j.resize(n_max);
  h.resize(n_max);

  sphericalBessel(r, n_max, j, y);

  for (n = 0; n < n_max; n++) {
    h1[n] = std::complex<double> (j[n], y[n]);
    h2[n] = std::complex<double> (j[n], -y[n]);
  }
}

// Calculate an - equation (5)
std::complex<double> calc_an(int n, double XL, std::complex<double> Ha, std::complex<double> mL,
	                         std::complex<double> PsiXL, std::complex<double> ZetaXL,
	                         std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1) {

  std::complex<double> Num = (Ha/mL + n/XL)*PsiXL - PsiXLM1;
  std::complex<double> Denom = (Ha/mL + n/XL)*ZetaXL - ZetaXLM1;

  return Num/Denom;
}

// Calculate bn - equation (6)
std::complex<double> calc_bn(int n, double XL, std::complex<double> Hb, std::complex<double> mL,
	                         std::complex<double> PsiXL, std::complex<double> ZetaXL,
	                         std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1) {

  std::complex<double> Num = (mL*Hb + n/XL)*PsiXL - PsiXLM1;
  std::complex<double> Denom = (mL*Hb + n/XL)*ZetaXL - ZetaXLM1;

  return Num/Denom;
}

// Calculates S1 - equation (25a)
std::complex<double> calc_S1(int n, std::complex<double> an, std::complex<double> bn,
		                     double Pi, double Tau) {

  return double(n + n + 1)*(Pi*an + Tau*bn)/double(n*n + n);
}

// Calculates S2 - equation (25b) (it's the same as (25a), just switches Pi and Tau)
std::complex<double> calc_S2(int n, std::complex<double> an, std::complex<double> bn,
				             double Pi, double Tau) {

  return calc_S1(n, an, bn, Tau, Pi);
}


//**********************************************************************************//
// This function calculates the Riccati-Bessel functions (Psi and Zeta) for a       //
// real argument (x).                                                               //
// Equations (20a) - (21b)                                                          //
//                                                                                  //
// Input parameters:                                                                //
//   x: Real argument to evaluate Psi and Zeta                                      //
//   n_max: Maximum number of terms to calculate Psi and Zeta                       //
//                                                                                  //
// Output parameters:                                                               //
//   Psi, Zeta: Riccati-Bessel functions                                            //
//**********************************************************************************//
void calcPsiZeta(double x, int n_max,
		         std::vector<std::complex<double> > D1,
		         std::vector<std::complex<double> > D3,
		         std::vector<std::complex<double> > &Psi,
		         std::vector<std::complex<double> > &Zeta) {

  int n;

  //Upward recurrence for Psi and Zeta - equations (20a) - (21b)
  Psi[0] = std::complex<double>(sin(x), 0);
  Zeta[0] = std::complex<double>(sin(x), -cos(x));
  for (n = 1; n <= n_max; n++) {
    Psi[n] = Psi[n - 1]*(n/x - D1[n - 1]);
    Zeta[n] = Zeta[n - 1]*(n/x - D3[n - 1]);
  }
}

//**********************************************************************************//
// This function calculates the logarithmic derivatives of the Riccati-Bessel       //
// functions (D1 and D3) for a complex argument (z).                                //
// Equations (16a), (16b) and (18a) - (18d)                                         //
//                                                                                  //
// Input parameters:                                                                //
//   z: Complex argument to evaluate D1 and D3                                      //
//   n_max: Maximum number of terms to calculate D1 and D3                          //
//                                                                                  //
// Output parameters:                                                               //
//   D1, D3: Logarithmic derivatives of the Riccati-Bessel functions                //
//**********************************************************************************//
void calcD1D3(std::complex<double> z, int n_max,
		      std::vector<std::complex<double> > &D1,
		      std::vector<std::complex<double> > &D3) {

  int n;
  std::vector<std::complex<double> > PsiZeta;
  PsiZeta.resize(n_max + 1);

  // Downward recurrence for D1 - equations (16a) and (16b)
  D1[n_max] = std::complex<double>(0.0, 0.0);
  for (n = n_max; n > 0; n--) {
    D1[n - 1] = double(n)/z - 1.0/(D1[n] + double(n)/z);
  }

  // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
  PsiZeta[0] = 0.5*(1.0 - std::complex<double>(cos(2.0*z.real()), sin(2.0*z.real()))*exp(-2.0*z.imag()));
  D3[0] = std::complex<double>(0.0, 1.0);
  for (n = 1; n <= n_max; n++) {
    PsiZeta[n] = PsiZeta[n - 1]*(double(n)/z - D1[n - 1])*(double(n)/z- D3[n - 1]);
    D3[n] = D1[n] + std::complex<double>(0.0, 1.0)/PsiZeta[n];
  }
}

//**********************************************************************************//
// This function calculates Pi and Tau for all values of Theta.                     //
// Equations (26a) - (26c)                                                          //
//                                                                                  //
// Input parameters:                                                                //
//   n_max: Maximum number of terms to calculate Pi and Tau                         //
//   nTheta: Number of scattering angles                                            //
//   Theta: Array containing all the scattering angles where the scattering         //
//          amplitudes will be calculated                                           //
//                                                                                  //
// Output parameters:                                                               //
//   Pi, Tau: Angular functions Pi and Tau, as defined in equations (26a) - (26c)   //
//**********************************************************************************//
void calcPiTau(int n_max, int nTheta, std::vector<double> Theta,
		       std::vector< std::vector<double> > &Pi,
		       std::vector< std::vector<double> > &Tau) {

  int n, t;
  for (n = 0; n < n_max; n++) {
    //****************************************************//
    // Equations (26a) - (26c)                            //
    //****************************************************//
    for (t = 0; t < nTheta; t++) {
      if (n == 0) {
        // Initialize Pi and Tau
        Pi[n][t] = 1.0;
        Tau[n][t] = (n + 1)*cos(Theta[t]);
      } else {
        // Calculate the actual values
        Pi[n][t] = ((n == 1) ? ((n + n + 1)*cos(Theta[t])*Pi[n - 1][t]/n)
                             : (((n + n + 1)*cos(Theta[t])*Pi[n - 1][t] - (n + 1)*Pi[n - 2][t])/n));
        Tau[n][t] = (n + 1)*cos(Theta[t])*Pi[n][t] - (n + 2)*Pi[n - 1][t];
      }
    }
  }
}

//**********************************************************************************//
// This function calculates the scattering coefficients required to calculate       //
// both the near- and far-field parameters.                                         //
//                                                                                  //
// Input parameters:                                                                //
//   L: Number of layers                                                            //
//   pl: Index of PEC layer. If there is none just send -1                          //
//   x: Array containing the size parameters of the layers [0..L-1]                 //
//   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
//   n_max: Maximum number of multipolar expansion terms to be used for the         //
//          calculations. Only used if you know what you are doing, otherwise set   //
//          this parameter to -1 and the function will calculate it.                //
//                                                                                  //
// Output parameters:                                                               //
//   an, bn: Complex scattering amplitudes                                          //
//                                                                                  //
// Return value:                                                                    //
//   Number of multipolar expansion terms used for the calculations                 //
//**********************************************************************************//
int ScattCoeffs(int L, int pl, std::vector<double> x, std::vector<std::complex<double> > m, int n_max,
		        std::vector<std::complex<double> > &an, std::vector<std::complex<double> > &bn) {
  //************************************************************************//
  // Calculate the index of the first layer. It can be either 0 (default)   //
  // or the index of the outermost PEC layer. In the latter case all layers //
  // below the PEC are discarded.                                           //
  //************************************************************************//

  int fl = (pl > 0) ? pl : 0;

  if (n_max <= 0) {
    n_max = Nmax(L, fl, pl, x, m);
  }

  std::complex<double> z1, z2;
  std::complex<double> Num, Denom;
  std::complex<double> G1, G2;
  std::complex<double> Temp;

  int n, l;

  //**************************************************************************//
  // Note that since Fri, Nov 14, 2014 all arrays start from 0 (zero), which  //
  // means that index = layer number - 1 or index = n - 1. The only exception //
  // are the arrays for representing D1, D3 and Q because they need a value   //
  // for the index 0 (zero), hence it is important to consider this shift     //
  // between different arrays. The change was done to optimize memory usage.  //
  //**************************************************************************//

  // Allocate memory to the arrays
  std::vector<std::vector<std::complex<double> > > D1_mlxl, D1_mlxlM1;
  D1_mlxl.resize(L);
  D1_mlxlM1.resize(L);

  std::vector<std::vector<std::complex<double> > > D3_mlxl, D3_mlxlM1;
  D3_mlxl.resize(L);
  D3_mlxlM1.resize(L);

  std::vector<std::vector<std::complex<double> > > Q;
  Q.resize(L);

  std::vector<std::vector<std::complex<double> > > Ha, Hb;
  Ha.resize(L);
  Hb.resize(L);

  for (l = 0; l < L; l++) {
    D1_mlxl[l].resize(n_max + 1);
    D1_mlxlM1[l].resize(n_max + 1);

    D3_mlxl[l].resize(n_max + 1);
    D3_mlxlM1[l].resize(n_max + 1);

    Q[l].resize(n_max + 1);

    Ha[l].resize(n_max);
    Hb[l].resize(n_max);
  }

  an.resize(n_max);
  bn.resize(n_max);

  std::vector<std::complex<double> > D1XL, D3XL;
  D1XL.resize(n_max + 1);
  D3XL.resize(n_max + 1);


  std::vector<std::complex<double> > PsiXL, ZetaXL;
  PsiXL.resize(n_max + 1);
  ZetaXL.resize(n_max + 1);

  //*************************************************//
  // Calculate D1 and D3 for z1 in the first layer   //
  //*************************************************//
  if (fl == pl) {  // PEC layer
    for (n = 0; n <= n_max; n++) {
      D1_mlxl[fl][n] = std::complex<double>(0.0, -1.0);
      D3_mlxl[fl][n] = std::complex<double>(0.0, 1.0);
    }
  } else { // Regular layer
    z1 = x[fl]* m[fl];

    // Calculate D1 and D3
    calcD1D3(z1, n_max, D1_mlxl[fl], D3_mlxl[fl]);
  }

  //******************************************************************//
  // Calculate Ha and Hb in the first layer - equations (7a) and (8a) //
  //******************************************************************//
  for (n = 0; n < n_max; n++) {
    Ha[fl][n] = D1_mlxl[fl][n + 1];
    Hb[fl][n] = D1_mlxl[fl][n + 1];
  }

  //*****************************************************//
  // Iteration from the second layer to the last one (L) //
  //*****************************************************//
  for (l = fl + 1; l < L; l++) {
    //************************************************************//
    //Calculate D1 and D3 for z1 and z2 in the layers fl+1..L     //
    //************************************************************//
    z1 = x[l]*m[l];
    z2 = x[l - 1]*m[l];

    //Calculate D1 and D3 for z1
    calcD1D3(z1, n_max, D1_mlxl[l], D3_mlxl[l]);

    //Calculate D1 and D3 for z2
    calcD1D3(z2, n_max, D1_mlxlM1[l], D3_mlxlM1[l]);

    //*********************************************//
    //Calculate Q, Ha and Hb in the layers fl+1..L //
    //*********************************************//

    // Upward recurrence for Q - equations (19a) and (19b)
    Num = exp(-2.0*(z1.imag() - z2.imag()))*std::complex<double>(cos(-2.0*z2.real()) - exp(-2.0*z2.imag()), sin(-2.0*z2.real()));
    Denom = std::complex<double>(cos(-2.0*z1.real()) - exp(-2.0*z1.imag()), sin(-2.0*z1.real()));
    Q[l][0] = Num/Denom;

    for (n = 1; n <= n_max; n++) {
      Num = (z1*D1_mlxl[l][n] + double(n))*(double(n) - z1*D3_mlxl[l][n - 1]);
      Denom = (z2*D1_mlxlM1[l][n] + double(n))*(double(n) - z2*D3_mlxlM1[l][n - 1]);

      Q[l][n] = (((x[l - 1]*x[l - 1])/(x[l]*x[l])* Q[l][n - 1])*Num)/Denom;
    }

    // Upward recurrence for Ha and Hb - equations (7b), (8b) and (12) - (15)
    for (n = 1; n <= n_max; n++) {
      //Ha
      if ((l - 1) == pl) { // The layer below the current one is a PEC layer
        G1 = -D1_mlxlM1[l][n];
        G2 = -D3_mlxlM1[l][n];
      } else {
        G1 = (m[l]*Ha[l - 1][n - 1]) - (m[l - 1]*D1_mlxlM1[l][n]);
        G2 = (m[l]*Ha[l - 1][n - 1]) - (m[l - 1]*D3_mlxlM1[l][n]);
      }

      Temp = Q[l][n]*G1;

      Num = (G2*D1_mlxl[l][n]) - (Temp*D3_mlxl[l][n]);
      Denom = G2 - Temp;

      Ha[l][n - 1] = Num/Denom;

      //Hb
      if ((l - 1) == pl) { // The layer below the current one is a PEC layer
        G1 = Hb[l - 1][n - 1];
        G2 = Hb[l - 1][n - 1];
      } else {
        G1 = (m[l - 1]*Hb[l - 1][n - 1]) - (m[l]*D1_mlxlM1[l][n]);
        G2 = (m[l - 1]*Hb[l - 1][n - 1]) - (m[l]*D3_mlxlM1[l][n]);
      }

      Temp = Q[l][n]*G1;

      Num = (G2*D1_mlxl[l][n]) - (Temp* D3_mlxl[l][n]);
      Denom = (G2- Temp);

      Hb[l][n - 1] = (Num/ Denom);
    }
  }

  //**************************************//
  //Calculate D1, D3, Psi and Zeta for XL //
  //**************************************//

  // Calculate D1XL and D3XL
  calcD1D3(x[L - 1], n_max, D1XL, D3XL);

  // Calculate PsiXL and ZetaXL
  calcPsiZeta(x[L - 1], n_max, D1XL, D3XL, PsiXL, ZetaXL);

  //*********************************************************************//
  // Finally, we calculate the scattering coefficients (an and bn) and   //
  // the angular functions (Pi and Tau). Note that for these arrays the  //
  // first layer is 0 (zero), in future versions all arrays will follow  //
  // this convention to save memory. (13 Nov, 2014)                      //
  //*********************************************************************//
  for (n = 0; n < n_max; n++) {
    //********************************************************************//
    //Expressions for calculating an and bn coefficients are not valid if //
    //there is only one PEC layer (ie, for a simple PEC sphere).          //
    //********************************************************************//
    if (pl < (L - 1)) {
      an[n] = calc_an(n + 1, x[L - 1], Ha[L - 1][n], m[L - 1], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
      bn[n] = calc_bn(n + 1, x[L - 1], Hb[L - 1][n], m[L - 1], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
    } else {
      an[n] = calc_an(n + 1, x[L - 1], std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
      bn[n] = PsiXL[n + 1]/ZetaXL[n + 1];
    }
  }

  return n_max;
}

//**********************************************************************************//
// This function calculates the actual scattering parameters and amplitudes         //
//                                                                                  //
// Input parameters:                                                                //
//   L: Number of layers                                                            //
//   pl: Index of PEC layer. If there is none just send -1                          //
//   x: Array containing the size parameters of the layers [0..L-1]                 //
//   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
//   nTheta: Number of scattering angles                                            //
//   Theta: Array containing all the scattering angles where the scattering         //
//          amplitudes will be calculated                                           //
//   n_max: Maximum number of multipolar expansion terms to be used for the         //
//          calculations. Only used if you know what you are doing, otherwise set   //
//          this parameter to -1 and the function will calculate it                 //
//                                                                                  //
// Output parameters:                                                               //
//   Qext: Efficiency factor for extinction                                         //
//   Qsca: Efficiency factor for scattering                                         //
//   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
//   Qbk: Efficiency factor for backscattering                                      //
//   Qpr: Efficiency factor for the radiation pressure                              //
//   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
//   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
//   S1, S2: Complex scattering amplitudes                                          //
//                                                                                  //
// Return value:                                                                    //
//   Number of multipolar expansion terms used for the calculations                 //
//**********************************************************************************//

int nMie(int L, int pl, std::vector<double> x, std::vector<std::complex<double> > m,
         int nTheta, std::vector<double> Theta, int n_max,
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
		 std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2)  {

  int i, n, t;
  std::vector<std::complex<double> > an, bn;
  std::complex<double> Qbktmp;

  // Calculate scattering coefficients
  n_max = ScattCoeffs(L, pl, x, m, n_max, an, bn);

  std::vector< std::vector<double> > Pi;
  Pi.resize(n_max);
  std::vector< std::vector<double> > Tau;
  Tau.resize(n_max);
  for (n = 0; n < n_max; n++) {
    Pi[n].resize(nTheta);
    Tau[n].resize(nTheta);
  }

  calcPiTau(n_max, nTheta, Theta, Pi, Tau);

  double x2 = x[L - 1]*x[L - 1];

  // Initialize the scattering parameters
  *Qext = 0;
  *Qsca = 0;
  *Qabs = 0;
  *Qbk = 0;
  Qbktmp = std::complex<double>(0.0, 0.0);
  *Qpr = 0;
  *g = 0;
  *Albedo = 0;

  // Initialize the scattering amplitudes
  for (t = 0; t < nTheta; t++) {
    S1[t] = std::complex<double>(0.0, 0.0);
    S2[t] = std::complex<double>(0.0, 0.0);
  }

  // By using downward recurrence we avoid loss of precision due to float rounding errors
  // See: https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
  //      http://en.wikipedia.org/wiki/Loss_of_significance
  for (i = n_max - 2; i >= 0; i--) {
    n = i + 1;
    // Equation (27)
    *Qext += (n + n + 1)*(an[i].real() + bn[i].real());
    // Equation (28)
    *Qsca += (n + n + 1)*(an[i].real()*an[i].real() + an[i].imag()*an[i].imag() + bn[i].real()*bn[i].real() + bn[i].imag()*bn[i].imag());
    // Equation (29)
    // We must check carefully this equation. If we remove the typecast to double then the result changes. Which is the correct one??? Ovidio (2014/12/10)
    *Qpr += ((n*(n + 2)/(n + 1))*((an[i]*std::conj(an[n]) + bn[i]*std::conj(bn[n])).real()) + ((double)(n + n + 1)/(n*(n + 1)))*(an[i]*std::conj(bn[i])).real());
    // Equation (33)
    Qbktmp = Qbktmp + (double)(n + n + 1)*(1 - 2*(n % 2))*(an[i]- bn[i]);

    //****************************************************//
    // Calculate the scattering amplitudes (S1 and S2)    //
    // Equations (25a) - (25b)                            //
    //****************************************************//
    for (t = 0; t < nTheta; t++) {
      S1[t] += calc_S1(n, an[i], bn[i], Pi[i][t], Tau[i][t]);
      S2[t] += calc_S2(n, an[i], bn[i], Pi[i][t], Tau[i][t]);
    }
  }

  *Qext = 2*(*Qext)/x2;                                 // Equation (27)
  *Qsca = 2*(*Qsca)/x2;                                 // Equation (28)
  *Qpr = *Qext - 4*(*Qpr)/x2;                           // Equation (29)

  *Qabs = *Qext - *Qsca;                                // Equation (30)
  *Albedo = *Qsca / *Qext;                              // Equation (31)
  *g = (*Qext - *Qpr) / *Qsca;                          // Equation (32)

  *Qbk = (Qbktmp.real()*Qbktmp.real() + Qbktmp.imag()*Qbktmp.imag())/x2;    // Equation (33)

  return n_max;
}

//**********************************************************************************//
// This function is just a wrapper to call the full 'nMie' function with fewer      //
// parameters, it is here mainly for compatibility with older versions of the       //
// program. Also, you can use it if you neither have a PEC layer nor want to define //
// any limit for the maximum number of terms.                                       //
//                                                                                  //
// Input parameters:                                                                //
//   L: Number of layers                                                            //
//   x: Array containing the size parameters of the layers [0..L-1]                 //
//   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
//   nTheta: Number of scattering angles                                            //
//   Theta: Array containing all the scattering angles where the scattering         //
//          amplitudes will be calculated                                           //
//                                                                                  //
// Output parameters:                                                               //
//   Qext: Efficiency factor for extinction                                         //
//   Qsca: Efficiency factor for scattering                                         //
//   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
//   Qbk: Efficiency factor for backscattering                                      //
//   Qpr: Efficiency factor for the radiation pressure                              //
//   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
//   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
//   S1, S2: Complex scattering amplitudes                                          //
//                                                                                  //
// Return value:                                                                    //
//   Number of multipolar expansion terms used for the calculations                 //
//**********************************************************************************//

int nMie(int L, std::vector<double> x, std::vector<std::complex<double> > m,
         int nTheta, std::vector<double> Theta,
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
         std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2) {

  return nMie(L, -1, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
}


//**********************************************************************************//
// This function is just a wrapper to call the full 'nMie' function with fewer      //
// parameters, it is useful if you want to include a PEC layer but not a limit      //
// for the maximum number of terms.                                                 //
//                                                                                  //
// Input parameters:                                                                //
//   L: Number of layers                                                            //
//   pl: Index of PEC layer. If there is none just send -1                          //
//   x: Array containing the size parameters of the layers [0..L-1]                 //
//   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
//   nTheta: Number of scattering angles                                            //
//   Theta: Array containing all the scattering angles where the scattering         //
//          amplitudes will be calculated                                           //
//                                                                                  //
// Output parameters:                                                               //
//   Qext: Efficiency factor for extinction                                         //
//   Qsca: Efficiency factor for scattering                                         //
//   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
//   Qbk: Efficiency factor for backscattering                                      //
//   Qpr: Efficiency factor for the radiation pressure                              //
//   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
//   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
//   S1, S2: Complex scattering amplitudes                                          //
//                                                                                  //
// Return value:                                                                    //
//   Number of multipolar expansion terms used for the calculations                 //
//**********************************************************************************//

int nMie(int L, int pl, std::vector<double> x, std::vector<std::complex<double> > m,
         int nTheta, std::vector<double> Theta,
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
         std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2) {

  return nMie(L, pl, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
}

//**********************************************************************************//
// This function is just a wrapper to call the full 'nMie' function with fewer      //
// parameters, it is useful if you want to include a limit for the maximum number   //
// of terms but not a PEC layer.                                                    //
//                                                                                  //
// Input parameters:                                                                //
//   L: Number of layers                                                            //
//   x: Array containing the size parameters of the layers [0..L-1]                 //
//   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
//   nTheta: Number of scattering angles                                            //
//   Theta: Array containing all the scattering angles where the scattering         //
//          amplitudes will be calculated                                           //
//   n_max: Maximum number of multipolar expansion terms to be used for the         //
//          calculations. Only used if you know what you are doing, otherwise set   //
//          this parameter to -1 and the function will calculate it                 //
//                                                                                  //
// Output parameters:                                                               //
//   Qext: Efficiency factor for extinction                                         //
//   Qsca: Efficiency factor for scattering                                         //
//   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
//   Qbk: Efficiency factor for backscattering                                      //
//   Qpr: Efficiency factor for the radiation pressure                              //
//   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
//   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
//   S1, S2: Complex scattering amplitudes                                          //
//                                                                                  //
// Return value:                                                                    //
//   Number of multipolar expansion terms used for the calculations                 //
//**********************************************************************************//

int nMie(int L, std::vector<double> x, std::vector<std::complex<double> > m,
         int nTheta, std::vector<double> Theta, int n_max,
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo,
         std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2) {

  return nMie(L, -1, x, m, nTheta, Theta, n_max, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
}


//**********************************************************************************//
// This function calculates complex electric and magnetic field in the surroundings //
// and inside (TODO) the particle.                                                  //
//                                                                                  //
// Input parameters:                                                                //
//   L: Number of layers                                                            //
//   pl: Index of PEC layer. If there is none just send 0 (zero)                    //
//   x: Array containing the size parameters of the layers [0..L-1]                 //
//   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
//   n_max: Maximum number of multipolar expansion terms to be used for the         //
//          calculations. Only used if you know what you are doing, otherwise set   //
//          this parameter to 0 (zero) and the function will calculate it.          //
//   nCoords: Number of coordinate points                                           //
//   Coords: Array containing all coordinates where the complex electric and        //
//           magnetic fields will be calculated                                     //
//                                                                                  //
// Output parameters:                                                               //
//   E, H: Complex electric and magnetic field at the provided coordinates          //
//                                                                                  //
// Return value:                                                                    //
//   Number of multipolar expansion terms used for the calculations                 //
//**********************************************************************************//

int nField(int L, int pl, std::vector<double> x, std::vector<std::complex<double> > m, int n_max,
           int nCoords, std::vector<double> Xp, std::vector<double> Yp, std::vector<double> Zp,
		   std::vector<std::complex<double> > &E, std::vector<std::complex<double> >  &H)  {

  int i, n, c;
  std::vector<std::complex<double> > an, bn;

  // Calculate scattering coefficients
  n_max = ScattCoeffs(L, pl, x, m, n_max, an, bn);

  std::vector< std::vector<double> > Pi, Tau;
  Pi.resize(n_max);
  Tau.resize(n_max);
  for (n = 0; n < n_max; n++) {
    Pi[n].resize(nCoords);
    Tau[n].resize(nCoords);
  }

  std::vector<double> Rho, Phi, Theta;
  Rho.resize(nCoords);
  Phi.resize(nCoords);
  Theta.resize(nCoords);

  for (c = 0; c < nCoords; c++) {
    // Convert to spherical coordinates
    Rho = sqrt(Xp[c]*Xp[c] + Yp[c]*Yp[c] + Zp[c]*Zp[c]);
    if (Rho < 1e-3) {
      Rho = 1e-3;
    }
    Phi = acos(Xp[c]/sqrt(Xp[c]*Xp[c] + Yp[c]*Yp[c]));
    Theta = acos(Xp[c]/Rho[c]);
  }

  calcPiTau(n_max, nCoords, Theta, Pi, Tau);

  std::vector<double > j, y;
  std::vector<std::complex<double> > h1, h2;
  j.resize(n_max);
  y.resize(n_max);
  h1.resize(n_max);
  h2.resize(n_max);

  for (c = 0; c < nCoords; c++) {
    //*******************************************************//
    // external scattering field = incident + scattered      //
    // BH p.92 (4.37), 94 (4.45), 95 (4.50)                  //
    // assume: medium is non-absorbing; refim = 0; Uabs = 0  //
    //*******************************************************//

    // Calculate spherical Bessel and Hankel functions
    sphericalBessel(Rho, n_max, j, y);
    sphericalHankel(Rho, n_max, h1, h2);

    // Initialize the fields
    E[c] = std::complex<double>(0.0, 0.0);
    H[c] = std::complex<double>(0.0, 0.0);

    // Firstly the easiest case, we want the field outside the particle
    if (Rho >= x[L - 1]) {
      
    }
  }

  return n_max;
}

