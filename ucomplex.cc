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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ucomplex.h"

complex Cadd(complex z1, complex z2) {
  complex zr;
  zr.r = z1.r + z2.r;
  zr.i = z1.i + z2.i;
  return zr;
}

complex Csub(complex z1, complex z2) {
  complex zr;
  zr.r = z1.r - z2.r;
  zr.i = z1.i - z2.i;
  return zr;
}

complex Cmul(complex z1, complex z2) {
  complex zr;
  zr.r = (z1.r*z2.r) - (z1.i*z2.i);
  zr.i = (z1.r*z2.i) + (z1.i*z2.r);
  return zr;
}

complex RCmul(double r, complex z) {
  complex zr;
  zr.r = z.r*r;
  zr.i = z.i*r;
  return zr;
}

complex Cdiv(complex z1, complex z2) {
/* The following algorithm is used to properly handle
   denominator overflow:

              |  a + b(d/c)   b - a(d/c)
              |  ---------- + ----------- I     if |d| < |c|
   a + b I    |  c + d(d/c)   c + d(d/c)
   -------  = |
   c + d I    |  b + a(c/d)   -a + b(c/d)
              |  ---------- + ----------- I     if |d| >= |c|
              |  d + c(c/d)   d + c(c/d)
*/
  complex zr;
  double tmp, denom;

  if(fabs(z2.r) > fabs(z2.i)) {
    tmp = z2.i/z2.r;
    denom = z2.r + z2.i*tmp;
    zr.r = (z1.r + z1.i*tmp)/denom;
    zr.i = (z1.i - z1.r*tmp)/denom;
  } else {
    tmp = z2.r/z2.i;
    denom = z2.i + z2.r*tmp;
    zr.r = (z1.i + z1.r*tmp)/denom;
    zr.i = (-z1.r + z1.i*tmp)/denom;
  }
  return zr;
}    

complex Complex(double re, double im) {
  complex zr;
  zr.r = re;
  zr.i = im;
  return zr;
}

double Cabs(complex z) {
  double r;
  r = sqrt((z.r*z.r) + (z.i*z.i));
  return r;
}

double Carc(complex z) {
  double r;
  r = atan2(z.i, z.r);
  return r;
}

complex Conjg(complex z) {
  complex zr;
  zr.r = z.r;
  zr.i = -z.i;
  return zr;
}

complex Cinv(complex z) {
  complex zr;
  double denom;

  denom = (z.r*z.r) + (z.i*z.i);
  zr.r = z.r/denom;
  zr.i = -z.i/denom;
  return zr;
}

complex Cexp(complex z) {
// exp(a + ib) = exp(a)*exp(ib) = exp(a)*[cos(b) + i sin(b)]
  complex zr;
  double expz;

  expz = exp(z.r);
  zr.r = expz*cos(z.i);
  zr.i = expz*sin(z.i);
  return zr;
}

complex Clog(complex z) {
// ln( p exp(i0)) = ln(p) + i0 + 2kpi
  complex zr;
  zr.r = log(Cabs(z));
  zr.i = atan2(z.i, z.r);
  return zr;
}

complex Csqrt(complex z) {
  complex zr;
  double root, q;

  if((z.r != 0.0) ||(z.i != 0.0)) {
    root = sqrt(0.5*(fabs(z.r) + Cabs(z)));
    q = z.i/(2.0*root);
    if(z.r >= 0.0) {
      zr.r = root;
      zr.i = q;
    } else if(z.i < 0.0) {
             zr.r = -q;
             zr.i = -root;
           } else {
             zr.r =  q;
             zr.i =  root;
           }
    } else zr = z;
  return zr;
}

// complex trigonometric functions

// complex cosinus
complex Ccos(complex z) {
// cos(x+iy) = cos(x).cos(iy) - sin(x).sin(iy)
// cos(ix) = cosh(x) et sin(ix) = i.sinh(x)
  complex zr;
  zr.r = cos(z.r)*cosh(z.i);
  zr.i = -sin(z.r)*sinh(z.i);
  return zr;
}

// complex sinus
complex Csin(complex z) {
// sin(x+iy) = sin(x).cos(iy) + cos(x).sin(iy)
// cos(ix) = cosh(x) et sin(ix) = i.sinh(x)
  complex zr;
  zr.r = sin(z.r)*cosh(z.i);
  zr.i = cos(z.r)*sinh(z.i);
  return zr;
}

// tangent
complex Ctan(complex z) {
  return Cdiv(Csin(z), Ccos(z));
}

// Inverse trigonometric functions

// complex arc cosinus
complex Carc_cos(complex z) {
// arccos(z) = -i.arcch(z)
  complex temp = Carc_ch(z);
  return Complex(temp.i, -temp.r);
}

// complex arc sinus
complex Carc_sin(complex z) {
// arcsin(z) = -i.arcsh(i.z)
  complex temp = Carc_sh(Complex(-z.i, z.r));
  return Complex(temp.i, -temp.r);
}

// complex arc tangent
complex Carc_tan(complex z) {
// arctg(z) = -i.arcth(i.z)
  complex temp = Carc_th(Complex(-z.i, z.r));
  return Complex(temp.i, -temp.r);
}

// Hyberbolic complex functions

// complex hyberbolic cosinus
complex Cch(complex z) {
// cosh(x+iy) = cosh(x).cosh(iy) + sinh(x).sinh(iy)
// cosh(iy) = cos(y) et sinh(iy) = i.sin(y)
  complex zr;
  zr.r = cosh(z.r)*cos(z.i);
  zr.i = sinh(z.r)*sin(z.i);
  return zr;
}

// complex hyberbolic sinus
complex Csh(complex z) {
// sinh(x+iy) = sinh(x).cosh(iy) + cosh(x).sinh(iy)
// cosh(iy) = cos(y) et sinh(iy) = i.sin(y)
  complex zr;
  zr.r = sinh(z.r)*cos(z.i);
  zr.i = cosh(z.r)*sin(z.i);
  return zr;
}

// complex hyberbolic tangent
complex Cth(complex z) {
// th(x) = sinh(x)/cosh(x)
// cosh(x) > 1 qq x
  complex temp, zr;
  temp = Cch(z);
  zr = Csh(z);
  return Cdiv(zr, temp);
}

// Inverse complex hyperbolic functions

// complex hyperbolic arc cosinus
complex Carc_ch(complex z) {
//                        _________
// arcch(z) = -/+ ln(z + V z.z - 1 )
  complex zr;
  complex z2 = Cmul(z, z);
  zr = Clog(Cadd(z, Csqrt(Complex(z2.r - 1.0, z2.i))));
  zr.r = -zr.r;
  zr.i = -zr.i;
  return zr;
}

// complex hyperbolic arc sinus
complex Carc_sh(complex z) {
//                    ________
// arcsh(z) = ln(z + V 1 + z.z)
  complex z2 = Cmul(z, z);
  return Clog(Cadd(z, Csqrt(Complex(z2.r + 1.0, z2.i))));
}

// complex hyperbolic arc tangent
complex Carc_th(complex z) {
// arcth(z) = 1/2 ln((z + 1)/(1 - z))
  return RCmul(0.5, Clog(Cdiv(Complex(z.r + 1.0, z.i), Complex(1.0 - z.r, -z.i))));
}

