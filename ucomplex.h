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

typedef struct COMPLEX {double r,i;} complex;

// Complex functions
complex Cadd(complex z1, complex z2);  // addition
complex Csub(complex z1, complex z2);  // substraction
complex Cmul(complex z1, complex z2);  // multiplication
complex RCmul(double r, complex z);    // double*complex
complex Cdiv(complex z1, complex z2);  // division
complex Complex(double r, double i);   // convert to double

complex Conjg(complex z);    // conjuge
complex Cinv(complex z);     // inverse function 1/z

// Complex functions with double return values
double Cabs(complex z);      // module
double Carc(complex z);      // argument : z1 / z = p.e^ia

// Elementary functions
complex Cexp(complex z);    // exponential
complex Clog(complex z);    // natural logarithm
complex Csqrt(complex z);   // square root

// Complex trigonometric functions 
complex Ccos(complex z);    // cosinus
complex Csin(complex z);    // sinus
complex Ctan(complex z);    // tangent

// Inverse complex trigonometric functions
complex Carc_cos(complex z);   // arc cosinus
complex Carc_sin(complex z);   // arc sinus
complex Carc_tan(complex z);   // arc tangent

// Hyperbolic complex functions
complex Cch(complex z);     // hyperbolic cosinus
complex Csh(complex z);     // hyperbolic sinus
complex Cth(complex z);     // hyperbolic tangent

// Inverse hyperbolic complex functions
complex Carc_ch(complex z); // hyperbolic arc cosinus
complex Carc_sh(complex z); // hyperbolic arc sinus
complex Carc_th(complex z); // hyperbolic arc tangente


