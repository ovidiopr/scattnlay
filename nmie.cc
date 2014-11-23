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
#include "ucomplex.h"
#include "nmie.h"

// TODO should be replaced with std::round() or std::lround()
#define round(x) ((x) >= 0 ? (int)((x) + 0.5):(int)((x) - 0.5))

complex C_ZERO = {0.0, 0.0};
complex C_ONE = {1.0, 0.0};
complex C_I = {0.0, 1.0};

int compare(std::string operation,  complex *a, std::vector< std::complex<double> > b) {
  for (int i = 0; i < b.size(); ++i) {
    //if (i > 50) continue;
    double diff_r = std::abs((a[i].r - b[i].real())/a[i].r);
    double diff_i = std::abs((a[i].i - b[i].imag())/a[i].i);
    double epsilon= 1e-16;
    if (diff_r > epsilon ||diff_i > epsilon) {
      printf("\n*** WARNING!! Non-zero diff!!! ***\n");
      printf("Op: %s at i=%i, diff_r=%g, diff_i=%g, a_r=%g, a_i = %g\n",
    	     operation.c_str(), i, diff_r, diff_i, a[i].r, a[i].i);
    }
    double factor = 1.0;
    if ((diff_r > epsilon * factor && a[i].r/a[0].r > epsilon)
	|| (diff_i > epsilon*factor && a[i].i/a[0].i > epsilon)) {
      printf("\n******** ERROR!! Non-zero diff!!! ********\n");
      printf("Op: %s at i=%i, diff_r=%g, diff_i=%g, a_r=%g, a_i = %g\n",
	     operation.c_str(), i, diff_r, diff_i, a[i].r, a[i].i);
    }
  }
  return 0;
}

int firstLayer(int L, int pl) {
  if (pl >= 0) {
    return pl;
  } else {
    return 0;
  }
}

// Calculate Nstop - equation (17)
int Nstop(double xL) {
  int result;

  if (xL <= 8) {
    result = std::lround(xL + 4*pow(xL, 1/3) + 1);
  } else if (xL <= 4200) {
    result = std::lround(xL + 4.05*pow(xL, 1/3) + 2);
  } else {
    result = std::lround(xL + 4*pow(xL, 1/3) + 2);
  }

  return result;
}

int Nmax(int L, int fl, int pl, double x[], complex m[]) {
  int i, result, ri, riM1;

  result = Nstop(x[L - 1]);
  for (i = fl; i < L; i++) {
    if (i > pl) {
      ri = round(Cabs(RCmul(x[i], m[i])));
    } else {
      ri = 0;
    }
    if (result < ri) {
      result = ri;
    }

    if ((i > fl) && ((i - 1) > pl)) {
      riM1 = round(Cabs(RCmul(x[i - 1], m[i])));
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
int Nmax_std(int L, int fl, int pl, std::vector<double> x,
		   std::vector<std::complex<double> > m) {
  int i, result, ri, riM1;
  result = Nstop(x[L - 1]);
  for (i = fl; i < L; i++) {
    if (i > pl) {
      ri = std::lround(std::abs(x[i]*m[i]));
    } else {
      ri = 0;
    }
    if (result < ri) {
      result = ri;
    }

    if ((i > fl) && ((i - 1) > pl)) {
      riM1 = std::lround(std::abs(x[i - 1]* m[i]));
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
// This function calculates the spherical Bessel functions (jn and hn) for a given  //
// value of r.                                                                      //
//                                                                                  //
// Input parameters:                                                                //
//   r: Real argument to evaluate jn and hn                                         //
//   n_max: Maximum number of terms to calculate jn and hn                          //
//                                                                                  //
// Output parameters:                                                               //
//   jn, hn: Spherical Bessel functions (double)                                    //
//**********************************************************************************//
void sphericalBessel(double r, int n_max, double **j, double **h) {
  int n;

  (*j)[0] = sin(r)/r;
  (*j)[1] = sin(r)/r/r - cos(r)/r;
  (*h)[0] = -cos(r)/r;
  (*h)[1] = -cos(r)/r/r - sin(r)/r;
  for (n = 2; n < n_max; n++) {
    (*j)[n] = (n + n + 1)*(*j)[n - 1]/r - (*j)[n - 2];
    (*h)[n] = (n + n + 1)*(*h)[n - 1]/r - (*h)[n - 2];
  }
}

//**********************************************************************************//
// This function calculates the spherical Bessel functions (jn and hn) for a given  //
// value of z.                                                                      //
//                                                                                  //
// Input parameters:                                                                //
//   z: Real argument to evaluate jn and hn                                         //
//   n_max: Maximum number of terms to calculate jn and hn                          //
//                                                                                  //
// Output parameters:                                                               //
//   jn, hn: Spherical Bessel functions (complex)                                   //
//**********************************************************************************//
void CsphericalBessel(complex z, int n_max, complex **j, complex **h) {
  int n;

  (*j)[0] = Cdiv(Csin(z), z);
  (*j)[1] = Csub(Cdiv(Cdiv(Csin(z), z), z), Cdiv(Ccos(z), z));
  (*h)[0] = Csub(C_ZERO, Cdiv(Ccos(z), z));
  (*h)[1] = Csub(C_ZERO, Cadd(Cdiv(Cdiv(Ccos(z), z), z), Cdiv(Csin(z), z)));
  for (n = 2; n < n_max; n++) {
    (*j)[n] = Csub(RCmul(n + n + 1, Cdiv((*j)[n - 1], z)), (*j)[n - 2]);
    (*h)[n] = Csub(RCmul(n + n + 1, Cdiv((*h)[n - 1], z)), (*h)[n - 2]);
  }
}

// Calculate an - equation (5)
complex calc_an(int n, double XL, complex Ha, complex mL, complex PsiXL, complex ZetaXL, complex PsiXLM1, complex ZetaXLM1) {
  complex Num = Csub(Cmul(Cadd(Cdiv(Ha, mL), Complex(n/XL, 0)), PsiXL), PsiXLM1);
  complex Denom = Csub(Cmul(Cadd(Cdiv(Ha, mL), Complex(n/XL, 0)), ZetaXL), ZetaXLM1);

  return Cdiv(Num, Denom);
}
/////////////////////////////////////////////
std::complex<double>
calc_an_std(int n, double XL, std::complex<double> Ha, std::complex<double> mL,
	    std::complex<double> PsiXL, std::complex<double> ZetaXL,
	    std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1) {
  std::complex<double> Num = (( (Ha/mL) + std::complex<double>(n/XL, 0) ) * PsiXL) - PsiXLM1;
  std::complex<double> Denom = (( (Ha/mL) + std::complex<double>(n/XL, 0) ) * ZetaXL) - ZetaXLM1;
  return Num/Denom;
}

// Calculate bn - equation (6)
complex calc_bn(int n, double XL, complex Hb, complex mL, complex PsiXL, complex ZetaXL, complex PsiXLM1, complex ZetaXLM1) {
  complex Num = Csub(Cmul(Cadd(Cmul(Hb, mL), Complex(n/XL, 0)), PsiXL), PsiXLM1);
  complex Denom = Csub(Cmul(Cadd(Cmul(Hb, mL), Complex(n/XL, 0)), ZetaXL), ZetaXLM1);

  return Cdiv(Num, Denom);
}
// Calculate bn - equation (6)
std::complex<double>
calc_bn_std(int n, double XL, std::complex<double> Hb, std::complex<double> mL,
	    std::complex<double> PsiXL, std::complex<double> ZetaXL,
	    std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1) {
 std::complex<double> Num = (( (Hb*mL) + std::complex<double>(n/XL,0) ) * PsiXL) - PsiXLM1;
 std::complex<double> Denom = (( (Hb*mL) + std::complex<double>(n/XL,0) ) * ZetaXL) - ZetaXLM1;

  return Num/Denom;
}

// Calculates S1_n - equation (25a)
complex calc_S1_n(int n, complex an, complex bn, double Pin, double Taun) {
  return RCmul((double)(n + n + 1)/(double)(n*n + n), Cadd(RCmul(Pin, an), RCmul(Taun, bn)));
}
//////////////////////
std::complex<double> calc_S1_n_std(int n, std::complex<double> an, std::complex<double> bn,
		      double Pin, double Taun) {
  return (double)(n + n + 1)/(double)(n*n + n) * ((Pin*an) + (Taun*bn));
}

// Calculates S2_n - equation (25b) (it's the same as (25a), just switches Pin and Taun)
complex calc_S2_n(int n, complex an, complex bn, double Pin, double Taun) {
  return calc_S1_n(n, an, bn, Taun, Pin);
}
//////////////////////
std::complex<double> calc_S2_n_std(int n, std::complex<double> an, std::complex<double> bn,
				   double Pin, double Taun) {
  return calc_S1_n_std(n, an, bn, Taun, Pin);
}


//**********************************************************************************//
// This function calculates the Riccati-Bessel functions (Psi and Zeta) for a       //
// given value of z.                                                                //
//                                                                                  //
// Input parameters:                                                                //
//   z: Complex argument to evaluate Psi and Zeta                                   //
//   n_max: Maximum number of terms to calculate Psi and Zeta                       //
//                                                                                  //
// Output parameters:                                                               //
//   Psi, Zeta: Riccati-Bessel functions                                            //
//**********************************************************************************//
void calcPsiZeta(complex z, int n_max, complex *D1, complex *D3, complex **Psi, complex **Zeta) {
  int n;
  complex cn;

  //Upward recurrence for Psi and Zeta - equations (20a) - (21b)
  (*Psi)[0] = Complex(sin(z.r), 0);
  (*Zeta)[0] = Complex(sin(z.r), -cos(z.r));
  for (n = 1; n <= n_max; n++) {
    cn = Complex(n, 0);
    (*Psi)[n] = Cmul((*Psi)[n - 1], Csub(Cdiv(cn, z), D1[n - 1]));
    (*Zeta)[n] = Cmul((*Zeta)[n - 1], Csub(Cdiv(cn, z), D3[n - 1]));
  }
}
//**********************************************************************************//
void calcPsiZeta_std(std::complex<double> z, int n_max,
		     std::vector< std::complex<double> > D1,
		     std::vector< std::complex<double> > D3,
		     std::vector< std::complex<double> > &Psi,
		     std::vector< std::complex<double> > &Zeta) {
  int n;
  std::complex<double> cn;

  //Upward recurrence for Psi and Zeta - equations (20a) - (21b)
  Psi[0] = std::complex<double>(sin(z.real()), 0);
  Zeta[0] = std::complex<double>(sin(z.real()), -cos(z.real()));
  for (n = 1; n <= n_max; n++) {
    cn = std::complex<double>(n, 0);
    Psi[n] = Psi[n-1] * ( (cn/z) - D1[n-1] );
    Zeta[n] = Zeta[n-1] * ( (cn/z) - D3[n-1] );
  }
}

//**********************************************************************************//
// This function calculates the logarithmic derivatives of the Riccati-Bessel       //
// functions (D1 and D3) for a given value of z.                                    //
//                                                                                  //
// Input parameters:                                                                //
//   z: Complex argument to evaluate D1 and D3                                      //
//   n_max: Maximum number of terms to calculate D1 and D3                          //
//                                                                                  //
// Output parameters:                                                               //
//   D1, D3: Logarithmic derivatives of the Riccati-Bessel functions                //
//**********************************************************************************//
void calcD1D3(complex z, int n_max, complex **D1, complex **D3) {
  int n;
  complex cn;
  complex *PsiZeta = (complex *) malloc((n_max + 1)*sizeof(complex));

  // Downward recurrence for D1 - equations (16a) and (16b)
  (*D1)[n_max] = C_ZERO;
  for (n = n_max; n > 0; n--) {
    cn = Complex(n, 0);
    (*D1)[n - 1] = Csub(Cdiv(cn, z), Cdiv(C_ONE, Cadd((*D1)[n], Cdiv(cn, z))));
  }

  // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
  PsiZeta[0] = RCmul(0.5, Csub(C_ONE, Cmul(Complex(cos(2*z.r), sin(2*z.r)), Complex(exp(-2*z.i), 0))));
  (*D3)[0] = C_I;
  for (n = 1; n <= n_max; n++) {
    cn = Complex(n, 0);
    PsiZeta[n] = Cmul(PsiZeta[n - 1], Cmul(Csub(Cdiv(cn, z), (*D1)[n - 1]), Csub(Cdiv(cn, z), (*D3)[n - 1])));
    (*D3)[n] = Cadd((*D1)[n], Cdiv(C_I, PsiZeta[n]));
  }

  free(PsiZeta);
}
//**********************************************************************************//
void calcD1D3_std(std::complex<double> z, int n_max,
		  std::vector< std::complex<double> > &D1,
		  std::vector< std::complex<double> > &D3) {
  int n;
  std::complex<double> cn;
  //complex *PsiZeta = (complex *) malloc((n_max + 1)*sizeof(complex));
  std::vector< std::complex<double> > PsiZeta;
  PsiZeta.resize(n_max + 1);
  // Downward recurrence for D1 - equations (16a) and (16b)
  D1[n_max] = std::complex<double>(0.0, 0.0);
  for (n = n_max; n > 0; n--) {
    cn = std::complex<double>(n, 0.0); 
    D1[n - 1] = cn/z - 1.0/(D1[n] + (cn/z));
  }

  // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
  PsiZeta[0] = 0.5 * (1.0 - (std::complex<double>(cos(2*z.real()), sin(2*z.real()))
			     * std::complex<double>(exp(-2*z.imag()), 0)));
  D3[0] = std::complex<double>(0.0, 1.0);
  for (n = 1; n <= n_max; n++) {
    cn = std::complex<double>(n, 0.0); 
    PsiZeta[n] = PsiZeta[n-1] * (((cn/ z) - D1[n - 1]) * ((cn/ z)- D3[n - 1]));
    D3[n] = D1[n]+ (std::complex<double>(0.0, 1.0)/ PsiZeta[n]);
  }
  // free(PsiZeta);
}

//**********************************************************************************//
// This function calculates Pi and Tau for all values of Theta.                     //
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
void calcPiTau(int n_max, int nTheta, double Theta[], double ***Pi, double ***Tau) {
  int n, t;

  for (n = 0; n < n_max; n++) {
    //****************************************************//
    // Equations (26a) - (26c)                            //
    //****************************************************//
    for (t = 0; t < nTheta; t++) {
      if (n == 0) {
        // Initialize Pi and Tau
        (*Pi)[n][t] = 1.0;
        (*Tau)[n][t] = (n + 1)*cos(Theta[t]);
      } else {
        // Calculate the actual values
        (*Pi)[n][t] = ((n == 1) ?  ((n + n + 1)*cos(Theta[t])*(*Pi)[n - 1][t]/n)
                                : (((n + n + 1)*cos(Theta[t])*(*Pi)[n - 1][t] - (n + 1)*(*Pi)[n - 2][t])/n));
        (*Tau)[n][t] = (n + 1)*cos(Theta[t])*(*Pi)[n][t] - (n + 2)*(*Pi)[n - 1][t];
      }
    }
  }
}
//**********************************************************************************//
void calcPiTau_std(int n_max, int nTheta, std::vector<double> Theta,
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
        Pi[n][t] = ((n == 1) ?  ((n + n + 1)*cos(Theta[t])*Pi[n - 1][t]/n)
                                :
		    (((n + n + 1)*cos(Theta[t])*Pi[n - 1][t] - (n + 1)*Pi[n - 2][t])/n));
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
int ScattCoeff(int L, int pl, double x[], complex m[], int n_max, complex **an, complex **bn){
  //************************************************************************//
  // Calculate the index of the first layer. It can be either 0 (default)   //
  // or the index of the outermost PEC layer. In the latter case all layers //
  // below the PEC are discarded.                                           //
  //************************************************************************//
  int fl = firstLayer(L, pl);

  if (n_max <= 0) {
    n_max = Nmax(L, fl, pl, x, m);
  }
  
  complex z1, z2, cn;
  complex Num, Denom;
  complex G1, G2;
  complex Temp;
  double Tmp;

  int n, l, t;

  //**************************************************************************//
  // Note that since Fri, Nov 14, 2014 all arrays start from 0 (zero), which  //
  // means that index = layer number - 1 or index = n - 1. The only exception //
  // are the arrays for representing D1, D3 and Q because they need a value   //
  // for the index 0 (zero), hence it is important to consider this shift     //
  // between different arrays. The change was done to optimize memory usage.  //
  //**************************************************************************//

  // Allocate memory to the arrays
  complex **D1_mlxl = (complex **) malloc(L*sizeof(complex *));
  complex **D1_mlxlM1 = (complex **) malloc(L*sizeof(complex *));

  complex **D3_mlxl = (complex **) malloc(L*sizeof(complex *));
  complex **D3_mlxlM1 = (complex **) malloc(L*sizeof(complex *));

  complex **Q = (complex **) malloc(L*sizeof(complex *));

  complex **Ha = (complex **) malloc(L*sizeof(complex *));
  complex **Hb = (complex **) malloc(L*sizeof(complex *));

  for (l = 0; l < L; l++) {
    D1_mlxl[l] = (complex *) malloc((n_max + 1)*sizeof(complex));
    D1_mlxlM1[l] = (complex *) malloc((n_max + 1)*sizeof(complex));

    D3_mlxl[l] = (complex *) malloc((n_max + 1)*sizeof(complex));
    D3_mlxlM1[l] = (complex *) malloc((n_max + 1)*sizeof(complex));

    Q[l] = (complex *) malloc((n_max + 1)*sizeof(complex));

    Ha[l] = (complex *) malloc(n_max*sizeof(complex));
    Hb[l] = (complex *) malloc(n_max*sizeof(complex));
  }

  (*an) = (complex *) malloc(n_max*sizeof(complex));
  (*bn) = (complex *) malloc(n_max*sizeof(complex));

  complex *D1XL = (complex *) malloc((n_max + 1)*sizeof(complex));
  complex *D3XL = (complex *) malloc((n_max + 1)*sizeof(complex));

  complex *PsiXL = (complex *) malloc((n_max + 1)*sizeof(complex));
  complex *ZetaXL = (complex *) malloc((n_max + 1)*sizeof(complex));

  //*************************************************//
  // Calculate D1 and D3 for z1 in the first layer   //
  //*************************************************//
  if (fl == pl) {  // PEC layer
    for (n = 0; n <= n_max; n++) {
      D1_mlxl[fl][n] = Complex(0, -1);
      D3_mlxl[fl][n] = C_I;
    }
  } else { // Regular layer
    z1 = RCmul(x[fl], m[fl]);

    // Calculate D1 and D3
    calcD1D3(z1, n_max, &(D1_mlxl[fl]), &(D3_mlxl[fl]));
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
    z1 = RCmul(x[l], m[l]);
    z2 = RCmul(x[l - 1], m[l]);

    //Calculate D1 and D3 for z1
    calcD1D3(z1, n_max, &(D1_mlxl[l]), &(D3_mlxl[l]));

    //Calculate D1 and D3 for z2
    calcD1D3(z2, n_max, &(D1_mlxlM1[l]), &(D3_mlxlM1[l]));

    //*********************************************//
    //Calculate Q, Ha and Hb in the layers fl+1..L //
    //*********************************************//

    // Upward recurrence for Q - equations (19a) and (19b)
    Num = RCmul(exp(-2*(z1.i - z2.i)), Complex(cos(-2*z2.r) - exp(-2*z2.i), sin(-2*z2.r)));
    Denom = Complex(cos(-2*z1.r) - exp(-2*z1.i), sin(-2*z1.r));
    Q[l][0] = Cdiv(Num, Denom);

    for (n = 1; n <= n_max; n++) {
      cn = Complex(n, 0);
      Num = Cmul(Cadd(Cmul(z1, D1_mlxl[l][n]), cn), Csub(cn, Cmul(z1, D3_mlxl[l][n - 1])));
      Denom = Cmul(Cadd(Cmul(z2, D1_mlxlM1[l][n]), cn), Csub(cn, Cmul(z2, D3_mlxlM1[l][n - 1])));

      Q[l][n] = Cdiv(Cmul(RCmul((x[l - 1]*x[l - 1])/(x[l]*x[l]), Q[l][n - 1]), Num), Denom);
    }

    // Upward recurrence for Ha and Hb - equations (7b), (8b) and (12) - (15)
    for (n = 1; n <= n_max; n++) {
      //Ha
      if ((l - 1) == pl) { // The layer below the current one is a PEC layer
        G1 = RCmul(-1.0, D1_mlxlM1[l][n]);
        G2 = RCmul(-1.0, D3_mlxlM1[l][n]);
      } else {
        G1 = Csub(Cmul(m[l], Ha[l - 1][n - 1]), Cmul(m[l - 1], D1_mlxlM1[l][n]));
        G2 = Csub(Cmul(m[l], Ha[l - 1][n - 1]), Cmul(m[l - 1], D3_mlxlM1[l][n]));
      }

      Temp = Cmul(Q[l][n], G1);

      Num = Csub(Cmul(G2, D1_mlxl[l][n]), Cmul(Temp, D3_mlxl[l][n]));
      Denom = Csub(G2, Temp);

      Ha[l][n - 1] = Cdiv(Num, Denom);

      //Hb
      if ((l - 1) == pl) { // The layer below the current one is a PEC layer
        G1 = Hb[l - 1][n - 1];
        G2 = Hb[l - 1][n - 1];
      } else {
        G1 = Csub(Cmul(m[l - 1], Hb[l - 1][n - 1]), Cmul(m[l], D1_mlxlM1[l][n]));
        G2 = Csub(Cmul(m[l - 1], Hb[l - 1][n - 1]), Cmul(m[l], D3_mlxlM1[l][n]));
      }

      Temp = Cmul(Q[l][n], G1);

      Num = Csub(Cmul(G2, D1_mlxl[l][n]), Cmul(Temp, D3_mlxl[l][n]));
      Denom = Csub(G2, Temp);

      Hb[l][n - 1] = Cdiv(Num, Denom);
    }
  }

  //**************************************//
  //Calculate D1, D3, Psi and Zeta for XL //
  //**************************************//
  z1 = Complex(x[L - 1], 0);

  // Calculate D1XL and D3XL
  calcD1D3(z1, n_max, &D1XL, &D3XL);

  // Calculate PsiXL and ZetaXL
  calcPsiZeta(z1, n_max, D1XL, D3XL, &PsiXL, &ZetaXL);

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
      (*an)[n] = calc_an(n + 1, x[L - 1], Ha[L - 1][n], m[L - 1], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
      (*bn)[n] = calc_bn(n + 1, x[L - 1], Hb[L - 1][n], m[L - 1], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
    } else {
      (*an)[n] = calc_an(n + 1, x[L - 1], C_ZERO, C_ONE, PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
      (*bn)[n] = Cdiv(PsiXL[n + 1], ZetaXL[n + 1]);
    }
  }

  // Free the memory used for the arrays
  for (l = 0; l < L; l++) {
    free(D1_mlxl[l]);
    free(D1_mlxlM1[l]);

    free(D3_mlxl[l]);
    free(D3_mlxlM1[l]);

    free(Q[l]);

    free(Ha[l]);
    free(Hb[l]);
  }

  free(D1_mlxl);
  free(D1_mlxlM1);

  free(D3_mlxl);
  free(D3_mlxlM1);

  free(Q);

  free(Ha);
  free(Hb);

  free(D1XL);
  free(D3XL);

  free(PsiXL);
  free(ZetaXL);

  return n_max;
}
//**********************************************************************************//
//**********************************************************************************//
//**********************************************************************************//
int ScattCoeff_std(double x[], complex m[], 
		   int L, int pl, std::vector<double> x_std,
		   std::vector<std::complex<double> > m_std, int n_max,
		   std::vector< std::complex<double> > &an_std, 
		   std::vector< std::complex<double> > &bn_std){
  //************************************************************************//
  // Calculate the index of the first layer. It can be either 0 (default)   //
  // or the index of the outermost PEC layer. In the latter case all layers //
  // below the PEC are discarded.                                           //
  //************************************************************************//
  
  //TODO: Why?
  //  int fl = firstLayer(L, pl);
  // instead of
  int fl = (pl > 0) ? pl : 0;
  // fl - first layer, pl - pec layer.

  if (n_max <= 0) {
    int tmp_n_max = Nmax(L, fl, pl, x, m);
    n_max = Nmax_std(L, fl, pl, x_std, m_std);
    if (n_max != tmp_n_max) printf("n_max mismatch 2\n");
  }
  
  // complex z1, z2, cn;
  // complex Num, Denom;
  // complex G1, G2;
  // complex Temp;
  std::complex<double> z1, z2, cn;
  std::complex<double> Num, Denom;
  std::complex<double> G1, G2;
  std::complex<double> Temp;

  double Tmp;

  int n, l, t;

  //**************************************************************************//
  // Note that since Fri, Nov 14, 2014 all arrays start from 0 (zero), which  //
  // means that index = layer number - 1 or index = n - 1. The only exception //
  // are the arrays for representing D1, D3 and Q because they need a value   //
  // for the index 0 (zero), hence it is important to consider this shift     //
  // between different arrays. The change was done to optimize memory usage.  //
  //**************************************************************************//

  // Allocate memory to the arrays
  //complex **D1_mlxl = (complex **) malloc(L*sizeof(complex *));
  //complex **D1_mlxlM1 = (complex **) malloc(L*sizeof(complex *));
  std::vector< std::vector< std::complex<double> > > D1_mlxl;
  D1_mlxl.resize(L);
  std::vector< std::vector< std::complex<double> > >D1_mlxlM1;
  D1_mlxlM1.resize(L);

  // complex **D3_mlxl = (complex **) malloc(L*sizeof(complex *));
  // complex **D3_mlxlM1 = (complex **) malloc(L*sizeof(complex *));
  std::vector< std::vector< std::complex<double> > > D3_mlxl;
  D3_mlxl.resize(L);
  std::vector< std::vector< std::complex<double> > > D3_mlxlM1;
  D3_mlxlM1.resize(L);

  //complex **Q = (complex **) malloc(L*sizeof(complex *));
  std::vector< std::vector< std::complex<double> > > Q;
  Q.resize(L);

  //complex **Ha = (complex **) malloc(L*sizeof(complex *));
  std::vector< std::vector< std::complex<double> > > Ha;
  Ha.resize(L);
  //complex **Hb = (complex **) malloc(L*sizeof(complex *));
  std::vector< std::vector< std::complex<double> > > Hb;
  Hb.resize(L);

  for (l = 0; l < L; l++) {
    // D1_mlxl[l] = (complex *) malloc((n_max + 1)*sizeof(complex));
    // D1_mlxlM1[l] = (complex *) malloc((n_max + 1)*sizeof(complex));
    D1_mlxl[l].resize(n_max +1);
    D1_mlxlM1[l].resize(n_max +1);

    // D3_mlxl[l] = (complex *) malloc((n_max + 1)*sizeof(complex));
    // D3_mlxlM1[l] = (complex *) malloc((n_max + 1)*sizeof(complex));
    D3_mlxl[l].resize(n_max +1);
    D3_mlxlM1[l].resize(n_max +1);

    //Q[l] = (complex *) malloc((n_max + 1)*sizeof(complex));
    Q[l].resize(n_max + 1);

    // Ha[l] = (complex *) malloc(n_max*sizeof(complex));
    // Hb[l] = (complex *) malloc(n_max*sizeof(complex));
    Ha[l].resize(n_max);
    Hb[l].resize(n_max);
  }

  // (*an) = (complex *) malloc(n_max*sizeof(complex));
  // (*bn) = (complex *) malloc(n_max*sizeof(complex));
  an_std.resize(n_max);
  bn_std.resize(n_max);

  // complex *D1XL = (complex *) malloc((n_max + 1)*sizeof(complex));
  // complex *D3XL = (complex *) malloc((n_max + 1)*sizeof(complex));
  std::vector< std::complex<double> > D1XL;
  D1XL.resize(n_max+1);
  std::vector< std::complex<double> > D3XL;
  D3XL.resize(n_max+1);


  // complex *PsiXL = (complex *) malloc((n_max + 1)*sizeof(complex));
  // complex *ZetaXL = (complex *) malloc((n_max + 1)*sizeof(complex));
  std::vector< std::complex<double> > PsiXL;
  PsiXL.resize(n_max+1);
  std::vector< std::complex<double> > ZetaXL;
  ZetaXL.resize(n_max+1);

  //*************************************************//
  // Calculate D1 and D3 for z1 in the first layer   //
  //*************************************************//
  if (fl == pl) {  // PEC layer
    for (n = 0; n <= n_max; n++) {
      D1_mlxl[fl][n] = std::complex<double>(0, -1);
      D3_mlxl[fl][n] = std::complex<double>(0, 1);
    }
  } else { // Regular layer
    z1 = x_std[fl]* m_std[fl];

    // Calculate D1 and D3
    calcD1D3_std(z1, n_max, D1_mlxl[fl], D3_mlxl[fl]);
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
    z1 = x_std[l] * m_std[l];
    z2 = x_std[l - 1] * m_std[l];

    //Calculate D1 and D3 for z1
    calcD1D3_std(z1, n_max, D1_mlxl[l], D3_mlxl[l]);

    //Calculate D1 and D3 for z2
    calcD1D3_std(z2, n_max, D1_mlxlM1[l], D3_mlxlM1[l]);

    //*********************************************//
    //Calculate Q, Ha and Hb in the layers fl+1..L //
    //*********************************************//

    // Upward recurrence for Q - equations (19a) and (19b)
    //Num = RCmul(exp(-2*(z1.i - z2.i)), Complex(cos(-2*z2.r) - exp(-2*z2.i), sin(-2*z2.r)));
    Num = exp( -2.0 * ( z1.imag() - z2.imag() ) ) *
      std::complex<double>(cos(-2.0*z2.real()) - exp(-2.0*z2.imag()), sin(-2.0*z2.real()));
    Denom = std::complex<double>(cos(-2.0*z1.real()) - exp(-2.0*z1.imag()),
				 sin(-2.0*z1.real()));
    Q[l][0] = Num/Denom;

    for (n = 1; n <= n_max; n++) {
      cn = std::complex<double>(n, 0);
      Num = ( (z1*D1_mlxl[l][n]) + cn) * (cn - (z1*D3_mlxl[l][n - 1]) );
      Denom = ( (z2*D1_mlxlM1[l][n]) + cn) * (cn- (z2* D3_mlxlM1[l][n - 1]) );

      Q[l][n] = (((x[l - 1]*x[l - 1])/(x[l]*x[l])* Q[l][n - 1]) * Num)/Denom;
    }

    // Upward recurrence for Ha and Hb - equations (7b), (8b) and (12) - (15)
    for (n = 1; n <= n_max; n++) {
      //Ha
      if ((l - 1) == pl) { // The layer below the current one is a PEC layer
        G1 = -1.0 * D1_mlxlM1[l][n];
        G2 = -1.0 * D3_mlxlM1[l][n];
      } else {
        G1 = (m_std[l] * Ha[l - 1][n - 1]) - (m_std[l - 1] * D1_mlxlM1[l][n]);
        G2 = (m_std[l] * Ha[l - 1][n - 1]) - (m_std[l - 1] * D3_mlxlM1[l][n]);
      }

      Temp = Q[l][n] * G1;

      Num = (G2*D1_mlxl[l][n]) - (Temp*D3_mlxl[l][n]);
      Denom = G2 - Temp;

      Ha[l][n - 1] = Num / Denom;

      //Hb
      if ((l - 1) == pl) { // The layer below the current one is a PEC layer
        G1 = Hb[l - 1][n - 1];
        G2 = Hb[l - 1][n - 1];
      } else {
        G1 = (m_std[l - 1] * Hb[l - 1][n - 1]) - (m_std[l] * D1_mlxlM1[l][n]);
        G2 = (m_std[l - 1] * Hb[l - 1][n - 1]) - (m_std[l] * D3_mlxlM1[l][n]);
      }

      Temp = Q[l][n] * G1;

      Num = (G2*D1_mlxl[l][n]) - (Temp* D3_mlxl[l][n]);
      Denom = (G2- Temp);

      Hb[l][n - 1] = (Num/ Denom);
    }
  }

  //**************************************//
  //Calculate D1, D3, Psi and Zeta for XL //
  //**************************************//
  z1 = std::complex<double>(x_std[L - 1], 0);

  // Calculate D1XL and D3XL
  calcD1D3_std(z1, n_max, D1XL, D3XL);

  // Calculate PsiXL and ZetaXL
  calcPsiZeta_std(z1, n_max, D1XL, D3XL, PsiXL, ZetaXL);

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
      an_std[n] = calc_an_std(n + 1, x_std[L-1], Ha[L-1][n], m_std[L-1],
			  PsiXL[n + 1], ZetaXL[n+1], PsiXL[n], ZetaXL[n]);
      bn_std[n] = calc_bn_std(n + 1, x_std[L-1], Hb[L-1][n], m_std[L-1],
			  PsiXL[n+1], ZetaXL[n+1], PsiXL[n], ZetaXL[n]);
    } else {
      an_std[n] = calc_an_std(n+1, x_std[L-1], std::complex<double>(0.0,0.0),
			  std::complex<double>(1.0,0.0),
			  PsiXL[n+1], ZetaXL[n+1], PsiXL[n], ZetaXL[n]);
      bn_std[n] = PsiXL[n+1] / ZetaXL[n+1];
    }
  }

  return n_max;
}


//**********************************************************************************//
// This function is just a wrapper to call the function 'nMieScatt' with fewer      //
// parameters, it is here mainly for compatibility with older versions of the       //
// program. Also, you can use it if you neither have a PEC layer nor want to define //
// any limit for the maximum number of terms.                                       //
//**********************************************************************************//

int nMie(int L, double x[], complex m[], int nTheta, double Theta[], double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[]) {

  return nMieScatt(L, -1, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
}


int nMie_std(double x[], complex m[], double Theta[], complex S1[], complex S2[],
int L, std::vector<double> &x_std, std::vector<std::complex<double> > &m_std, int nTheta, std::vector<double> &Theta_std, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector< std::complex<double> > &S1_std, std::vector< std::complex<double> >  &S2_std) {
  return nMieScatt_std(x, m, Theta, S1, S2, L, -1, x_std, m_std, nTheta, Theta_std, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1_std, S2_std);
}


//**********************************************************************************//
// This function is just a wrapper to call the function 'nMieScatt' with fewer      //
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

int nMiePEC(int L, int pl, double x[], complex m[], int nTheta, double Theta[], double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[]) {

  return nMieScatt(L, pl, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
}

//**********************************************************************************//
// This function is just a wrapper to call the function 'nMieScatt' with fewer      //
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

int nMieMax(int L, double x[], complex m[], int nTheta, double Theta[], int n_max, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[]) {

  return nMieScatt(L, -1, x, m, nTheta, Theta, n_max, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
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

int nMieScatt(int L, int pl, double x[], complex m[], int nTheta, double Theta[], int n_max, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[]) {
  int i, n, t;
  double **Pi, **Tau;
  complex *an, *bn;
  complex Qbktmp;

  n_max = ScattCoeff(L, pl, x, m, n_max, &an, &bn);

  Pi = (double **) malloc(n_max*sizeof(double *));
  Tau = (double **) malloc(n_max*sizeof(double *));
  for (n = 0; n < n_max; n++) {
    Pi[n] = (double *) malloc(nTheta*sizeof(double));
    Tau[n] = (double *) malloc(nTheta*sizeof(double));
  }

  calcPiTau(n_max, nTheta, Theta, &Pi, &Tau);

  double x2 = x[L - 1]*x[L - 1];

  // Initialize the scattering parameters
  *Qext = 0;
  *Qsca = 0;
  *Qabs = 0;
  *Qbk = 0;
  Qbktmp = C_ZERO;
  *Qpr = 0;
  *g = 0;
  *Albedo = 0;

  // Initialize the scattering amplitudes
  for (t = 0; t < nTheta; t++) {
    S1[t] = C_ZERO;
    S2[t] = C_ZERO;
  }

  // By using downward recurrence we avoid loss of precision due to float rounding errors
  // See: https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
  //      http://en.wikipedia.org/wiki/Loss_of_significance
  for (i = n_max - 2; i >= 0; i--) {
    n = i + 1;
    // Equation (27)
    *Qext = *Qext + (double)(n + n + 1)*(an[i].r + bn[i].r);
    // Equation (28)
    *Qsca = *Qsca + (double)(n + n + 1)*(an[i].r*an[i].r + an[i].i*an[i].i + bn[i].r*bn[i].r + bn[i].i*bn[i].i);
    // Equation (29)
    *Qpr = *Qpr + ((n*(n + 2)/(n + 1))*((Cadd(Cmul(an[i], Conjg(an[n])), Cmul(bn[i], Conjg(bn[n])))).r) + ((double)(n + n + 1)/(n*(n + 1)))*(Cmul(an[i], Conjg(bn[i])).r));

    // Equation (33)
    Qbktmp = Cadd(Qbktmp, RCmul((double)((n + n + 1)*(1 - 2*(n % 2))), Csub(an[i], bn[i])));

    //****************************************************//
    // Calculate the scattering amplitudes (S1 and S2)    //
    // Equations (25a) - (25b)                            //
    //****************************************************//
    for (t = 0; t < nTheta; t++) {
      S1[t] = Cadd(S1[t], calc_S1_n(n, an[i], bn[i], Pi[i][t], Tau[i][t]));
      S2[t] = Cadd(S2[t], calc_S2_n(n, an[i], bn[i], Pi[i][t], Tau[i][t]));
    }
  }

  *Qext = 2*(*Qext)/x2;                                 // Equation (27)
  *Qsca = 2*(*Qsca)/x2;                                 // Equation (28)
  *Qpr = *Qext - 4*(*Qpr)/x2;                           // Equation (29)

  *Qabs = *Qext - *Qsca;                                // Equation (30)
  *Albedo = *Qsca / *Qext;                              // Equation (31)
  *g = (*Qext - *Qpr) / *Qsca;                          // Equation (32)

  *Qbk = (Qbktmp.r*Qbktmp.r + Qbktmp.i*Qbktmp.i)/x2;    // Equation (33)

  // Free the memory used for the arrays
  for (n = 0; n < n_max; n++) {
    free(Pi[n]);
    free(Tau[n]);
  }

  free(Pi);
  free(Tau);

  free(an);
  free(bn);

  return n_max;
}



int nMieScatt_std(double x[], complex m[], double Theta[], complex S1[], complex S2[],
		  int L, int pl,
		  std::vector<double> &x_std, std::vector<std::complex<double> > &m_std,
		  int nTheta, std::vector<double> &Theta_std,
		  int n_max, double *Qext, double *Qsca, double *Qabs, double *Qbk,
		  double *Qpr, double *g, double *Albedo,
		  std::vector< std::complex<double> > &S1_std,
		  std::vector< std::complex<double> >  &S2_std)  {
  int i, n, t;
  // double **Pi, **Tau;
  std::vector< std::vector<double> > Pi_std, Tau_std;
  complex *an, *bn;
  std::vector< std::complex<double> > an_std, bn_std;
  //complex Qbktmp;
  std::complex<double> Qbktmp;
  {
    int tmp_n_max = ScattCoeff(L, pl, x, m, n_max, &an, &bn);
    n_max = ScattCoeff_std(x, m, L, pl, x_std, m_std, n_max, an_std, bn_std);
    if (n_max != tmp_n_max) printf("n_max mismatch\n");
    compare("an vs an_std: ", an, an_std);
    compare("bn vs bn_std: ", an, an_std);
  }
  // Pi = (double **) malloc(n_max*sizeof(double *));  
  // Tau = (double **) malloc(n_max*sizeof(double *));
  std::vector< std::vector<double> > Pi;
  Pi.resize(n_max);
  std::vector< std::vector<double> > Tau;
  Tau.resize(n_max);
  for (n = 0; n < n_max; n++) {
    // Pi[n] = (double *) malloc(nTheta*sizeof(double));
    // Tau[n] = (double *) malloc(nTheta*sizeof(double));
    Pi[n].resize(nTheta);
    Tau[n].resize(nTheta);
  }

  calcPiTau_std(n_max, nTheta, Theta_std, Pi, Tau);

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
    S1_std[t] = std::complex<double>(0.0, 0.0);
    S2_std[t] = std::complex<double>(0.0, 0.0);
  }

  // By using downward recurrence we avoid loss of precision due to float rounding errors
  // See: https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
  //      http://en.wikipedia.org/wiki/Loss_of_significance
  for (i = n_max - 2; i >= 0; i--) {
    n = i + 1;
    // Equation (27)
    *Qext = *Qext + (double)(n + n + 1)*(an_std[i].real() + bn_std[i].real());
    // Equation (28)
    *Qsca = *Qsca + (double)(n + n + 1)*(an_std[i].real()*an_std[i].real() + an_std[i].imag()*an_std[i].imag() + bn_std[i].real()*bn_std[i].real() + bn_std[i].imag()*bn_std[i].imag());
    // Equation (29)
  // Qpr = *Qpr + ((n*(n + 2)/(n + 1))*((Cadd(Cmul(an[i], Conjg(an[n])), Cmul(bn[i], Conjg(bn[n])))).r) + ((double)(n + n + 1)/(n*(n + 1)))*(Cmul(an[i], Conjg(bn[i])).r));
    *Qpr = *Qpr +
      (
       (n*(n + 2)/(n + 1))
       *(
	 (
	  (an_std[i]*std::conj(an_std[n]))
	  + (bn_std[i]*std::conj(bn_std[n]))
	  ).real()
	 )
       + ((double)(n + n + 1)/(n*(n + 1)))
	  *(
	    (an_std[i] * std::conj(bn_std[i])
	     ).real()
	    )
       );
    // Equation (33)
    Qbktmp = Qbktmp + ((double)((n + n + 1)*(1 - 2*(n % 2))) * (an_std[i]- bn_std[i]));

    //****************************************************//
    // Calculate the scattering amplitudes (S1 and S2)    //
    // Equations (25a) - (25b)                            //
    //****************************************************//
    for (t = 0; t < nTheta; t++) {
      S1_std[t] = S1_std[t] + calc_S1_n_std(n, an_std[i], bn_std[i], Pi[i][t], Tau[i][t]);
      S2_std[t] = S2_std[t] + calc_S2_n_std(n, an_std[i], bn_std[i], Pi[i][t], Tau[i][t]);
    }
  }

  *Qext = 2*(*Qext)/x2;                                 // Equation (27)
  *Qsca = 2*(*Qsca)/x2;                                 // Equation (28)
  *Qpr = *Qext - 4*(*Qpr)/x2;                           // Equation (29)

  *Qabs = *Qext - *Qsca;                                // Equation (30)
  *Albedo = *Qsca / *Qext;                              // Equation (31)
  *g = (*Qext - *Qpr) / *Qsca;                          // Equation (32)

  *Qbk = (Qbktmp.real()*Qbktmp.real() + Qbktmp.imag()*Qbktmp.imag())/x2;    // Equation (33)

  // Free the memory used for the arrays
  // for (n = 0; n < n_max; n++) {
  //   free(Pi[n]);
  //   free(Tau[n]);
  // }

  // free(Pi);
  // free(Tau);

  free(an);
  free(bn);

  return n_max;
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

int nMieField(int L, int pl, double x[], complex m[], int n_max, int nCoords, double Xp[], double Yp[], double Zp[], complex E[], complex H[]){
  int i, n, c;
  double **Pi, **Tau;
  complex *an, *bn;
  double *Rho = (double *) malloc(nCoords*sizeof(double));
  double *Phi = (double *) malloc(nCoords*sizeof(double));
  double *Theta = (double *) malloc(nCoords*sizeof(double));

  for (c = 0; c < nCoords; c++) {
    Rho[c] = sqrt(Xp[c]*Xp[c] + Yp[c]*Yp[c] + Zp[c]*Zp[c]);
    if (Rho[c] < 1e-3) {
      Rho[c] = 1e-3;
    }
    Phi[c] = acos(Xp[c]/sqrt(Xp[c]*Xp[c] + Yp[c]*Yp[c]));
    Theta[c] = acos(Xp[c]/Rho[c]);
  }

  n_max = ScattCoeff(L, pl, x, m, n_max, &an, &bn);

  Pi = (double **) malloc(n_max*sizeof(double *));
  Tau = (double **) malloc(n_max*sizeof(double *));
  for (n = 0; n < n_max; n++) {
    Pi[n] = (double *) malloc(nCoords*sizeof(double));
    Tau[n] = (double *) malloc(nCoords*sizeof(double));
  }

  calcPiTau(n_max, nCoords, Theta, &Pi, &Tau);

  double x2 = x[L - 1]*x[L - 1];

  // Initialize the fields
  for (c = 0; c < nCoords; c++) {
    E[c] = C_ZERO;
    H[c] = C_ZERO;
  }

  //*******************************************************//
  // external scattering field = incident + scattered      //
  // BH p.92 (4.37), 94 (4.45), 95 (4.50)                  //
  // assume: medium is non-absorbing; refim = 0; Uabs = 0  //
  //*******************************************************//

  // Firstly the easiest case, we want the field outside the particle
  if (Rho[c] >= x[L - 1]) {
  }
//  for (i = 1; i < (n_max - 1); i++) {
//    n = i - 1;
/*    // Equation (27)
    *Qext = *Qext + (double)(n + n + 1)*(an[i].r + bn[i].r);
    // Equation (28)
    *Qsca = *Qsca + (double)(n + n + 1)*(an[i].r*an[i].r + an[i].i*an[i].i + bn[i].r*bn[i].r + bn[i].i*bn[i].i);
    // Equation (29)
    *Qpr = *Qpr + ((n*(n + 2)/(n + 1))*((Cadd(Cmul(an[i], Conjg(an[n])), Cmul(bn[i], Conjg(bn[n])))).r) + ((double)(n + n + 1)/(n*(n + 1)))*(Cmul(an[i], Conjg(bn[i])).r));

    // Equation (33)
    Qbktmp = Cadd(Qbktmp, RCmul((double)((n + n + 1)*(1 - 2*(n % 2))), Csub(an[i], bn[i])));
*/
    //****************************************************//
    // Calculate the scattering amplitudes (S1 and S2)    //
    // Equations (25a) - (25b)                            //
    //****************************************************//
/*    for (t = 0; t < nTheta; t++) {
      S1[t] = Cadd(S1[t], calc_S1_n(n, an[i], bn[i], Pi[i][t], Tau[i][t]));
      S2[t] = Cadd(S2[t], calc_S2_n(n, an[i], bn[i], Pi[i][t], Tau[i][t]));
    }*/
//  }

  // Free the memory used for the arrays
  for (n = 0; n < n_max; n++) {
    free(Pi[n]);
    free(Tau[n]);
  }

  free(Pi);
  free(Tau);

  free(an);
  free(bn);

  free(Rho);
  free(Phi);
  free(Theta);

  return n_max;
}

