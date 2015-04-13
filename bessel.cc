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
#include "bessel.h"
#include <algorithm>
#include <complex>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nmie {
  namespace bessel {
    // Implementation of spherical Bessel functions from
    // L.-W. Cai / Computer Physics Communications 182 (2011) 663–668
    std::vector< std::complex<double> > bessel_j(int nmax, std::complex<double> rho) {
      if (nmax < 0) throw std::invalid_argument("Bessel order should be >= 0 (nmie::bessel_j)\n");
      std::complex<double> c_zero(0.0,0.0);
      std::vector< std::complex<double> > j(nmax+1, c_zero);
      // Select normalization amplitude
      int norm_n = 0;
      std::complex<double> norm_value = std::sin(rho)/rho,
	tmp = std::sin(rho)/(rho*rho) - std::cos(rho)/rho;
      if (std::abs(tmp) > std::abs(norm_value)) {
	norm_value = tmp;
	norm_n = 1;
      }
      // calculate number of terms for backward recursion Cai eq(11)
      double theta = std::arg(rho), r = std::abs(rho);
      double M = (1.83 + 4.1* std::pow(std::sin(theta), 0.36)) *
	std::pow(r, (0.91 - 0.43*std::pow(std::sin(theta),0.33)))
	+ 9.0*(1.0-std::sqrt(std::sin(theta)));
      int terms = std::max(nmax+2, static_cast<int>(std::ceil(M))+2);
      printf("terms = %d\n", terms);
      // Seed for eq(2)
      std::complex<double> b_np1(0.0, 0.0), b_nm1;
      double eps = std::numeric_limits<double>::epsilon();
      std::complex<double> b_n(eps, 0.0);
      //recurence
      for (int n = terms*2; n > 0; --n) {
	b_nm1 = (2.0*static_cast<double>(n)+1.0)/rho*b_n - b_np1;
	if (n-1 < nmax+1) j[n-1] = b_nm1;
      }
      // normalize
      norm_value /= j[norm_n];
      for (auto& value : j) value *=norm_value;
      
      return j;
    }
    // Implementation of Bessel functions from Bohren and Huffman book, pp. 86-87, eq 4.11
    // Calculate all orders of function from 0 to nmax (included) for argument rho
    std::vector< std::complex<double> > bh_bessel_j(int nmax, std::complex<double> rho) {
      if (nmax < 0) throw std::invalid_argument("Bessel order should be >= 0 (nmie::bessel_j)\n");
      std::vector< std::complex<double> > j(nmax+1);
      j[0] = std::sin(rho)/rho;
      if (nmax == 0) return j;
      j[1] = std::sin(rho)/(rho*rho) - std::cos(rho)/rho;
      if (nmax == 1) return j;
      for (int i = 2; i < nmax+1; ++i) {
	int n = i - 1;
	j[n+1] = static_cast<double>(2*n+1)/rho*j[n] - j[n-1];
      }
      return j;
    }


// !*****************************************************************************80
//  
//  C++ port of fortran code
//  
// !! CSPHJY: spherical Bessel functions jn(z) and yn(z) for complex argument.
// !
// !  Discussion:
// !
// !    This procedure computes spherical Bessel functions jn(z) and yn(z)
// !    and their derivatives for a complex argument.
// !
// !  Licensing:
// !
// !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
// !    they give permission to incorporate this routine into a user program 
// !    provided that the copyright is acknowledged.
// !
// !  Modified:
// !
// !    01 August 2012
// !
// !  Author:
// !
// !    Shanjie Zhang, Jianming Jin
// !
// !  Reference:
// !
// !    Shanjie Zhang, Jianming Jin,
// !    Computation of Special Functions,
// !    Wiley, 1996,
// !    ISBN: 0-471-11963-6,
// !    LC: QA351.C45.
// ! 
// !  Parameters:
// !
// !    Input, integer ( kind = 4 ) N, the order of jn(z) and yn(z).
// !
// !    Input, complex ( kind = 8 ) Z, the argument.
// !
// !    Output, integer ( kind = 4 ) NM, the highest order computed.
// !
// !    Output, complex ( kind = 8 ) CSJ(0:N0, CDJ(0:N), CSY(0:N), CDY(0:N),
// !    the values of jn(z), jn'(z), yn(z), yn'(z).
// !
    void csphjy (int n, std::complex<double>z, int& nm,
		 std::vector< std::complex<double> >& csj,
		 std::vector< std::complex<double> >& cdj,
		 std::vector< std::complex<double> >& csy,
		 std::vector< std::complex<double> >& cdy ) {
      double a0;
      csj.resize(n+1);
      cdj.resize(n+1);
      csy.resize(n+1);
      cdy.resize(n+1);
      std::complex<double> cf, cf0, cf1, cs, csa, csb;
      int k, m;
      a0 = std::abs(z);
      nm = n;      
      if (a0 < 1.0e-60) {
	for (int k = 0; k < n+1; ++k) {
	  csj[k] = 0.0;
	  cdj[k] = 0.0;
	  csy[k] = -1.0e+300;
	  cdy[k] = 1.0e+300;
	}
	csj[0] = std::complex<double>( 1.0, 0.0);
	cdj[1] =  std::complex<double>( 0.333333333333333, 0.0);
	return;
      }
      csj[0] = std::sin ( z ) / z;
      csj[1] = ( csj[0] - std::cos ( z ) ) / z;
      
      if ( 2 <= n ) {
	csa = csj[0];
	csb = csj[1];
	m = msta1 ( a0, 200 );
	if ( m < n ) nm = m;
	else m = msta2 ( a0, n, 15 );
	cf0 = 0.0;
	cf1 = 1.0e-100;
	for (int k = m; k>=0; --k) {
	  cf = ( 2.0 * k + 3.0 ) * cf1 / z - cf0;
	  if ( k <= nm ) csj[k] = cf;
	  cf0 = cf1;
	  cf1 = cf;
	}
	if ( std::abs ( csa ) <= std::abs ( csb ) ) cs = csb / cf0;
	else  cs = csa / cf;
	for (int k = 0; k <= nm; ++k) {
	  csj[k] = cs * csj[k];
	}
      }
      cdj[0] = ( std::cos ( z ) - std::sin ( z ) / z ) / z;
      for (int k = 1; k <=nm; ++k) {
	cdj[k] = csj[k-1] - ( k + 1.0 ) * csj[k] / z;
      }
      csy[0] = - std::cos ( z ) / z;
      csy[1] = ( csy[0] - std::sin ( z ) ) / z;
      cdy[0] = ( std::sin ( z ) + std::cos ( z ) / z ) / z;
      cdy[1] = ( 2.0 * cdy[0] - std::cos ( z ) )  / z;
      for (int k = 2; k<=nm; ++k) {
	if ( std::abs ( csj[k-2] ) < std::abs ( csj[k-1] ) ) {
	  csy[k] = ( csj[k] * csy[k-1] - 1.0 / ( z * z ) ) / csj[k-1];
	} else {
	  csy[k] = ( csj[k] * csy[k-2] - ( 2.0 * k - 1.0 ) / std::pow(z,3) )
	    / csj[k-2];
	}
      }
      for (int k = 2; k<=nm; ++k) {
	cdy[k] = csy[k-1] - ( k + 1.0 ) * csy[k] / z;
      }
      
      return;
    }
      // function msta2 ( x, n, mp )
      
      // !*****************************************************************************80
      // !
      // !! MSTA2 determines a backward recurrence starting point for Jn(x).
      // !
      // !  Discussion:
      // !
      // !    This procedure determines the starting point for a backward
      // !    recurrence such that all Jn(x) has MP significant digits.
      // !
      // !  Licensing:
      // !
      // !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
      // !    they give permission to incorporate this routine into a user program 
      // !    provided that the copyright is acknowledged.
      // !
      // !  Modified:
      // !
      // !    08 July 2012
      // !
      // !  Author:
      // !
      // !    Shanjie Zhang, Jianming Jin
      // !
      // !  Reference:
      // !
      // !    Shanjie Zhang, Jianming Jin,
      // !    Computation of Special Functions,
      // !    Wiley, 1996,
      // !    ISBN: 0-471-11963-6,
      // !    LC: QA351.C45.
      // !
      // !  Parameters:
      // !
      // !    Input, real ( kind = 8 ) X, the argument of Jn(x).
      // !
      // !    Input, integer ( kind = 4 ) N, the order of Jn(x).
      // !
      // !    Input, integer ( kind = 4 ) MP, the number of significant digits.
      // !
      // !    Output, integer ( kind = 4 ) MSTA2, the starting point.
      // !
    int msta2 ( double x, int n, int mp ) {
      double a0, ejn, f, f0, f1, hmp;
      int it, n0, n1, nn;
      double obj;
      a0 = std::abs ( x );
      hmp = 0.5 * mp;
      ejn = envj ( n, a0 );
      if ( ejn <= hmp ) {
	obj = mp;
	n0 = static_cast<int> ( 1.1 * a0 );
      } else {
	obj = hmp + ejn;
	n0 = n;
      }
      f0 = envj ( n0, a0 ) - obj;
      n1 = n0 + 5;
      f1 = envj ( n1, a0 ) - obj;
      for (int it = 1; it < 21; ++it) {
	nn = n1 - ( n1 - n0 ) / ( 1.0 - f0 / f1 );
	f = envj ( nn, a0 ) - obj;
	if ( std::abs ( nn - n1 ) < 1 ) break;
	n0 = n1;
	f0 = f1;
	n1 = nn;
	f1 = f;
      }
      return  nn + 10;
    }



// !*****************************************************************************80
// !
// !! MSTA1 determines a backward recurrence starting point for Jn(x).
// !
// !  Discussion:
// !
// !    This procedure determines the starting point for backward  
// !    recurrence such that the magnitude of    
// !    Jn(x) at that point is about 10^(-MP).
// !
// !  Licensing:
// !
// !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
// !    they give permission to incorporate this routine into a user program 
// !    provided that the copyright is acknowledged.
// !
// !  Modified:
// !
// !    08 July 2012
// !
// !  Author:
// !
// !    Shanjie Zhang, Jianming Jin
// !
// !  Reference:
// !
// !    Shanjie Zhang, Jianming Jin,
// !    Computation of Special Functions,
// !    Wiley, 1996,
// !    ISBN: 0-471-11963-6,
// !    LC: QA351.C45.
// !
// !  Parameters:
// !
// !    Input, real ( kind = 8 ) X, the argument.
// !
// !    Input, integer ( kind = 4 ) MP, the negative logarithm of the 
// !    desired magnitude.
// !
// !    Output, integer ( kind = 4 ) MSTA1, the starting point.
// !
    int  msta1 ( double x, int mp ) {
      double a0, f, f0, f1;
      int it, n0, n1, nn;
      a0 = std::abs ( x );
      n0 = static_cast<int> (1.1 * a0 ) + 1;
      f0 = envj ( n0, a0 ) - mp;
      n1 = n0 + 5;
      f1 = envj ( n1, a0 ) - mp;
      for (int it = 1; it <= 20; ++it) { 
	nn = n1 - ( n1 - n0 ) / ( 1.0 - f0 / f1 );
	f = envj ( nn, a0 ) - mp;
	if ( abs ( nn - n1 ) < 1 ) break;
	n0 = n1;
	f0 = f1;
	n1 = nn;
	f1 = f;
      }
      return nn;
    }
      // function envj ( n, x )

      // !*****************************************************************************80
      // !
      // !! ENVJ is a utility function used by MSTA1 and MSTA2.
      // !
      // !  Licensing:
      // !
      // !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
      // !    they give permission to incorporate this routine into a user program 
      // !    provided that the copyright is acknowledged.
      // !
      // !  Modified:
      // !
      // !    14 March 2012
      // !
      // !  Author:
      // !
      // !    Shanjie Zhang, Jianming Jin
      // !
      // !  Reference:
      // !
      // !    Shanjie Zhang, Jianming Jin,
      // !    Computation of Special Functions,
      // !    Wiley, 1996,
      // !    ISBN: 0-471-11963-6,
      // !    LC: QA351.C45.
      // !
      // !  Parameters:
      // !
      // !    Input, integer ( kind = 4 ) N, ?
      // !
      // !    Input, real ( kind = 8 ) X, ?
      // !
      // !    Output, real ( kind = 8 ) ENVJ, ?
      // !
    double envj (int n, double x ) {
      return  0.5 * std::log10(6.28 * n) - n * std::log10(1.36 * x / static_cast<double>(n) );
    }
  }  // end of namespace bessel
}  // end of namespace nmie
