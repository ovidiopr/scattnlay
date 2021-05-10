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
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "../../bessel.h"
#include "../../nmie.h"

int nmax_ = -1;
//Temporary variables
std::vector<std::complex<double> > PsiZeta_;

void calcD1D3(const std::complex<double> z, std::vector<std::complex<double> > &D1, std::vector<std::complex<double> > &D3) {
    D1[nmax_] = std::complex<double>(0.0, 0.0);
    const std::complex<double> zinv = std::complex<double>(1.0, 0.0)/z;
    for (int n = nmax_; n > 0; n--) {
      D1[n - 1] = static_cast<double>(n)*zinv - 1.0/(D1[n] + static_cast<double>(n)*zinv);
    }
    if (std::abs(D1[0]) > 100000.0)
      throw std::invalid_argument("Unstable D1! Please, try to change input parameters!\n");

    PsiZeta_.resize(nmax_ + 1);

    // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
    PsiZeta_[0] = 0.5*(1.0 - std::complex<double>(std::cos(2.0*z.real()), std::sin(2.0*z.real()))
                      *std::exp(-2.0*z.imag()));
    D3[0] = std::complex<double>(0.0, 1.0);
    for (int n = 1; n <= nmax_; n++) {
      PsiZeta_[n] = PsiZeta_[n - 1]*(static_cast<double>(n)*zinv - D1[n - 1])
                                   *(static_cast<double>(n)*zinv - D3[n - 1]);
      D3[n] = D1[n] + std::complex<double>(0.0, 1.0)/PsiZeta_[n];
    }
  }

  void calcPsiZeta(std::complex<double> z, std::vector<std::complex<double> > &Psi, std::vector<std::complex<double> > &Zeta) {
    std::complex<double> c_i(0.0, 1.0);
    std::vector<std::complex<double> > D1(nmax_ + 1), D3(nmax_ + 1);
    calcD1D3(z, D1, D3);
    Psi.resize(nmax_+1);
    Zeta.resize(nmax_+1);

    Psi[0] = std::sin(z);
    Zeta[0] = std::sin(z) - c_i*std::cos(z);
    for (int n = 1; n <= nmax_; n++) {
      Psi[n]  =  Psi[n - 1]*(static_cast<double>(n)/z - D1[n - 1]);
      Zeta[n] = Zeta[n - 1]*(static_cast<double>(n)/z - D3[n - 1]);
    }
  }


int main(int argc, char *argv[]) {

  int n_small = 1, n_big = 15;
  int n = n_big+2, n_print;
  std::complex<double> z_small(0.0012, 0.0034);
  std::complex<double> z_medium(1.0, 2.0);
  std::complex<double> z_big(18.0, 17.0);
  std::complex<double> z_big_real(81.0, 0.0);
  std::complex<double> z_big_imag(0.0, 13.0);
  std::complex<double> z;

  std::vector<std::complex<double> > D1, D3;

  int nm;
  // csj and cdj - complex argument spherical bessel first kind and derivative
  // csy and cdy - complex argument spherical bessel second kind and derivative
  std::vector< std::complex<double> > csj, cdj, csy, cdy, csh;
  std::vector<std::complex<double> > Psi, Zeta, dPsi, dZeta;

  auto c_i = std::complex<double>(0.0,1.0);

  printf("===== Small ========\n");
  z = z_small;
  printf("z = %+10.5er%+10.5ei\n", z.real(), z.imag());
  nmie::bessel::csphjy (n, z, nm, csj, cdj,  csy, cdy );
  csh.resize(csj.size());
  for (unsigned int i = 0; i < csj.size(); ++i)
    csh[i] = csj[i]+c_i*csy[i];
  nmie::bessel::calcPsi(n, z, Psi, dPsi);
  nmie::bessel::calcZeta(n, z, Zeta, dZeta);

  // n_print = n_small;
  // printf("csj[%i] = %+10.5er%+10.5ei;  csy[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csj[n_print].real(), csj[n_print].imag(), n_print,
  // 	 csy[n_print].real(), csy[n_print].imag());
  // printf("WA csj =         e   r        e   i;  WA csy =         e   r        e   i\n");
  // n_print = n_big;
  // printf("csj[%i] = %+10.5er%+10.5ei;  csy[%i] = %+10.5er%+10.5ei\n", n_big,
  // 	 csj[n_print].real(), csj[n_print].imag(), n_print,
  // 	 csy[n_print].real(), csy[n_print].imag());
  // printf(" WA csj =         e   r        e   i;   WA csy =         e   r        e   i\n");

  // n_print = n_small;
  // printf("csh[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csh[n_print].real(), csh[n_print].imag());
  // printf("WA csj =         e   r        e   i;\n");
  // n_print = n_big;
  // printf("csh[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csh[n_print].real(), csh[n_print].imag());
  // printf("WA csj  =         e   r        e   i;\n");

  n_print = n_small;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf("WA Psi =         e   r        e   i;  WA Zeta =         e   r        e   i\n");
  n_print = n_big;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf(" WA Psi =         e   r        e   i;   WA Zeta =         e   r        e   i\n");

  n_print = n_small;
  printf("dPsi[%i] = %+10.5er%+10.5ei;  dZeta[%i] = %+10.5er%+10.5ei\n", n_print,
	 dPsi[n_print].real(), dPsi[n_print].imag(), n_print,
	 dZeta[n_print].real(), dZeta[n_print].imag());
  printf("WA dPsi = 0.800005e-03r+0.226667e-02i;  WA dZeta = +4.8284 e+04r-5.98822e+04i\n");
  n_print = n_big;
  printf("dPsi[%i] = %+10.5er%+10.5ei;  dZeta[%i] = %+10.5er%+10.5ei\n", n_big,
	 dPsi[n_print].real(), dPsi[n_print].imag(), n_print,
	 dZeta[n_print].real(), dZeta[n_print].imag());
  printf(" WA dPsi = +1.75388e-53r-6.94429e-54i;   WA dZeta = +8.58553e+55r+7.47397e+55i\n");

  nmax_ = nm;
  printf("----- Scattnlay (nmax_=%d)-----\n",nmax_);
  calcPsiZeta(z, Psi, Zeta);

  n_print = n_small;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf("WA Psi =         e   r        e   i;  WA Zeta =         e   r        e   i\n");
  n_print = n_big;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf(" WA Psi =         e   r        e   i;   WA Zeta =         e   r        e   i\n");




  printf("===== Medium ========\n");
  z = z_medium;
  printf("z = %+10.5er%+10.5ei\n", z.real(), z.imag());
  nmie::bessel::csphjy (n, z, nm, csj, cdj,  csy, cdy );
  csh.resize(csj.size());
  for (unsigned int i = 0; i < csj.size(); ++i)
    csh[i] = csj[i]+c_i*csy[i];
  nmie::bessel::calcPsi(n, z, Psi, dPsi);
  nmie::bessel::calcZeta(n, z, Zeta, dZeta);

  // n_print = n_small;
  // printf("csj[%i] = %+10.5er%+10.5ei;  csy[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csj[n_print].real(), csj[n_print].imag(), n_print,
  // 	 csy[n_print].real(), csy[n_print].imag());
  // printf("WA csj =         e   r        e   i;  WA csy =         e   r        e   i\n");
  // n_print = n_big;
  // printf("csj[%i] = %+10.5er%+10.5ei;  csy[%i] = %+10.5er%+10.5ei\n", n_big,
  // 	 csj[n_print].real(), csj[n_print].imag(), n_print,
  // 	 csy[n_print].real(), csy[n_print].imag());
  // printf(" WA csj =         e   r        e   i;   WA csy =         e   r        e   i\n");

  // n_print = n_small;
  // printf("csh[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csh[n_print].real(), csh[n_print].imag());
  // printf("WA csj =         e   r        e   i;\n");
  // n_print = n_big;
  // printf("csh[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csh[n_print].real(), csh[n_print].imag());
  // printf("WA csj  =         e   r        e   i;\n");

  n_print = n_small;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf("WA Psi =         e   r        e   i;  WA Zeta =         e   r        e   i\n");
  n_print = n_big;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf(" WA Psi =         e   r        e   i;   WA Zeta =         e   r        e   i\n");

  n_print = n_small;
  printf("dPsi[%i] = %+10.5er%+10.5ei;  dZeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 dPsi[n_print].real(), dPsi[n_print].imag(), n_print,
  	 dZeta[n_print].real(), dZeta[n_print].imag());
  printf("WA dPsi = 2.41792e   r+1.27781 e   i;  WA dZeta = +1.99423e-01r-7.01483e-02i\n");
  n_print = n_big;
  printf("dPsi[%i] = %+10.5er%+10.5ei;  dZeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 dPsi[n_print].real(), dPsi[n_print].imag(), n_print,
  	 dZeta[n_print].real(), dZeta[n_print].imag());
  printf(" WA dPsi = -1.03313e-11r-1.13267e-11i;   WA dZeta = -2.11402e+11r+8.34528e+10i\n");

  nmax_ = nm;
  printf("----- Scattnlay (nmax_=%d)-----\n",nmax_);
  calcPsiZeta(z, Psi, Zeta);

  n_print = n_small;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf("WA Psi =         e   r        e   i;  WA Zeta =         e   r        e   i\n");
  n_print = n_big;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf(" WA Psi =         e   r        e   i;   WA Zeta =         e   r        e   i\n");

  printf("===== Big ========\n");
  z = z_big;
  printf("z = %+10.5er%+10.5ei\n", z.real(), z.imag());
  nmie::bessel::csphjy (n, z, nm, csj, cdj,  csy, cdy );
  csh.resize(csj.size());
  for (unsigned int i = 0; i < csj.size(); ++i)
    csh[i] = csj[i]+c_i*csy[i];
  nmie::bessel::calcPsi(n, z, Psi, dPsi);
  nmie::bessel::calcZeta(n, z, Zeta, dZeta);

  // n_print = n_small;
  // printf("csj[%i] = %+10.5er%+10.5ei;  csy[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csj[n_print].real(), csj[n_print].imag(), n_print,
  // 	 csy[n_print].real(), csy[n_print].imag());
  // printf("WA csj =         e   r        e   i;  WA csy =         e   r        e   i\n");
  // n_print = n_big;
  // printf("csj[%i] = %+10.5er%+10.5ei;  csy[%i] = %+10.5er%+10.5ei\n", n_big,
  // 	 csj[n_print].real(), csj[n_print].imag(), n_print,
  // 	 csy[n_print].real(), csy[n_print].imag());
  // printf(" WA csj =         e   r        e   i;   WA csy =         e   r        e   i\n");

  // n_print = n_small;
  // printf("csh[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csh[n_print].real(), csh[n_print].imag());
  // printf("WA csj =         e   r        e   i;\n");
  // n_print = n_big;
  // printf("csh[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csh[n_print].real(), csh[n_print].imag());
  // printf("WA csj  =         e   r        e   i;\n");

  n_print = n_small;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf("WA Psi =         e   r        e   i;  WA Zeta =         e   r        e   i\n");
  n_print = n_big;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf(" WA Psi =         e   r        e   i;   WA Zeta =         e   r        e   i\n");

  n_print = n_small;
  printf("dPsi[%i] = %+10.5er%+10.5ei;  dZeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 dPsi[n_print].real(), dPsi[n_print].imag(), n_print,
  	 dZeta[n_print].real(), dZeta[n_print].imag());
  printf("WA dPsi =         e   r        e   i;  WA dZeta = -3.11025e-08r-2.90558e-08i\n");
  n_print = n_big;
  printf("dPsi[%i] = %+10.5er%+10.5ei;  dZeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 dPsi[n_print].real(), dPsi[n_print].imag(), n_print,
  	 dZeta[n_print].real(), dZeta[n_print].imag());
  printf(" WA dPsi = -2.47740e+05r+3.05899e+05i;   WA dZeta = -6.13497e-07r-1.14967e-06i\n");

  nmax_ = 100;
  printf("----- Scattnlay (nmax_=%d)-----\n",nmax_);
  calcPsiZeta(z, Psi, Zeta);

  n_print = n_small;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf("WA Psi =         e   r        e   i;  WA Zeta =         e   r        e   i\n");
  n_print = n_big;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf(" WA Psi =         e   r        e   i;   WA Zeta =         e   r        e   i\n");


  printf("===== Big real========\n");
  z = z_big_real;
  printf("z = %+10.5er%+10.5ei\n", z.real(), z.imag());
  nmie::bessel::csphjy (n, z, nm, csj, cdj,  csy, cdy );
  csh.resize(csj.size());
  for (unsigned int i = 0; i < csj.size(); ++i)
    csh[i] = csj[i]+c_i*csy[i];
  nmie::bessel::calcPsi(n, z, Psi, dPsi);
  nmie::bessel::calcZeta(n, z, Zeta, dZeta);

  // n_print = n_small;
  // printf("csj[%i] = %+10.5er%+10.5ei;  csy[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csj[n_print].real(), csj[n_print].imag(), n_print,
  // 	 csy[n_print].real(), csy[n_print].imag());
  // printf("WA csj =         e   r        e   i;  WA csy =         e   r        e   i\n");
  // n_print = n_big;
  // printf("csj[%i] = %+10.5er%+10.5ei;  csy[%i] = %+10.5er%+10.5ei\n", n_big,
  // 	 csj[n_print].real(), csj[n_print].imag(), n_print,
  // 	 csy[n_print].real(), csy[n_print].imag());
  // printf(" WA csj =         e   r        e   i;   WA csy =         e   r        e   i\n");

  // n_print = n_small;
  // printf("csh[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csh[n_print].real(), csh[n_print].imag());
  // printf("WA csj =         e   r        e   i;\n");
  // n_print = n_big;
  // printf("csh[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csh[n_print].real(), csh[n_print].imag());
  // printf("WA csj  =         e   r        e   i;\n");

  n_print = n_small;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  //
  //SphericalHankelH1(1,81)*81
  printf("WA Psi =-0.784462377 e        e   i;  WA Zeta = -0.784462377r+0.620299278i\n");
  n_print = n_big;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf(" WA Psi = 0.699946236744       e   i;   WA Zeta = 0.69994623  r+0.727239996i\n");

  // n_print = n_small;
  // printf("dPsi[%i] = %+10.5er%+10.5ei;  dZeta[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 dPsi[n_print].real(), dPsi[n_print].imag(), n_print,
  // 	 dZeta[n_print].real(), dZeta[n_print].imag());
  // printf("WA dPsi =         e   r        e   i;  WA dZeta = -6.20203e-01r-7.84344e-01i\n");
  // n_print = n_big;
  // printf("dPsi[%i] = %+10.5er%+10.5ei;  dZeta[%i] = %+10.5er%+10.5ei\n", n_big,
  // 	 dPsi[n_print].real(), dPsi[n_print].imag(), n_print,
  // 	 dZeta[n_print].real(), dZeta[n_print].imag());
  // printf(" WA dPsi =         e   r        e   i;   WA dZeta = -7.13982e-01r+6.86858e-01i\n");

  nmax_ = 100 ;
  printf("----- Scattnlay (nmax_=%d)-----\n",nmax_);
  calcPsiZeta(z, Psi, Zeta);

  n_print = n_small;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf("WA Psi =         e   r        e   i;  WA Zeta =         e   r        e   i\n");
  n_print = n_big;
  printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  	 Zeta[n_print].real(), Zeta[n_print].imag());
  printf(" WA Psi =         e   r        e   i;   WA Zeta =         e   r        e   i\n");



  printf("===== Big imag========\n");
  z = z_big_imag;
  printf("z = %+10.5er%+10.5ei\n", z.real(), z.imag());
  nmie::bessel::csphjy (n, z, nm, csj, cdj,  csy, cdy );
  csh.resize(csj.size());
  for (unsigned int i = 0; i < csj.size(); ++i)
    csh[i] = csj[i]+c_i*csy[i];
  nmie::bessel::calcPsi(n, z, Psi, dPsi);
  nmie::bessel::calcZeta(n, z, Zeta, dZeta);

  // n_print = n_small;
  // printf("csj[%i] = %+10.5er%+10.5ei;  csy[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csj[n_print].real(), csj[n_print].imag(), n_print,
  // 	 csy[n_print].real(), csy[n_print].imag());
  // printf("WA csj =         e   r        e   i;  WA csy =         e   r        e   i\n");
  // n_print = n_big;
  // printf("csj[%i] = %+10.5er%+10.5ei;  csy[%i] = %+10.5er%+10.5ei\n", n_big,
  // 	 csj[n_print].real(), csj[n_print].imag(), n_print,
  // 	 csy[n_print].real(), csy[n_print].imag());
  // printf(" WA csj =         e   r        e   i;   WA csy =         e   r        e   i\n");

  // n_print = n_small;
  // printf("csh[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csh[n_print].real(), csh[n_print].imag());
  // printf("WA csj =         e   r        e   i;\n");
  // n_print = n_big;
  // printf("csh[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 csh[n_print].real(), csh[n_print].imag());
  // printf("WA csj  =         e   r        e   i;\n");

  // n_print = n_small;
  // printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_print,
  // 	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  // 	 Zeta[n_print].real(), Zeta[n_print].imag());
  // printf("WA Psi =         e   r        e   i;  WA Zeta =         e   r        e   i\n");
  // n_print = n_big;
  // printf("Psi[%i] = %+10.5er%+10.5ei;  Zeta[%i] = %+10.5er%+10.5ei\n", n_big,
  // 	 Psi[n_print].real(), Psi[n_print].imag(), n_print,
  // 	 Zeta[n_print].real(), Zeta[n_print].imag());
  // printf(" WA Psi =         e   r        e   i;   WA Zeta =         e   r        e   i\n");

  n_print = n_small;
  printf("dPsi[%i] = %+10.5er%+10.5ei;  dZeta[%i] = %+10.5er%+10.5ei\n", n_print,
  	 dPsi[n_print].real(), dPsi[n_print].imag(), n_print,
  	 dZeta[n_print].real(), dZeta[n_print].imag());
  printf("WA dPsi =         e   r        e   i;  WA dZeta = -2.47233e-22r-2.44758e-06i\n");
  n_print = n_big;
  printf("dPsi[%i] = %+10.5er%+10.5ei;  dZeta[%i] = %+10.5er%+10.5ei\n", n_big,
  	 dPsi[n_print].real(), dPsi[n_print].imag(), n_print,
  	 dZeta[n_print].real(), dZeta[n_print].imag());
  printf(" WA dPsi =         e   r        e   i;   WA dZeta =  8.39449e-19r+1.28767e-02i\n");

}

