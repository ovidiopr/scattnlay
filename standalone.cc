//**********************************************************************************//
//    Copyright (C) 2009-2013  Ovidio Pena <ovidio@bytesfall.com>                   //
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

#include <algorithm>
#include <complex>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "nmie.h"

const double PI=3.14159265358979323846;

//***********************************************************************************//
// This is the main function of 'scattnlay', here we read the parameters as          //
// arguments passed to the program which should be executed with the following       //
// syntaxis:                                                                         //
// ./scattnlay -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] [-t ti tf nt] [-c comment]  //
//                                                                                   //
// When all the parameters were correctly passed we setup the integer L (the         //
// number of layers) and the arrays x and m, containing the size parameters and      //
// refractive indexes of the layers, respectively and call the function nMie.        //
// If the calculation is successful the results are printed with the following       //
// format:                                                                           //
//                                                                                   //
//    * If no comment was passed:                                                    //
//        'Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo'                                    //
//                                                                                   //
//    * If a comment was passed:                                                     //
//        'comment, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo'                           //
//***********************************************************************************//
int main(int argc, char *argv[]) {
  try {
    char comment[200];
    int has_comment = 0;
    int i, l, L;
    std::vector<double> x, Theta;
    std::vector<std::complex<double> > m, S1, S2;
    double Qext, Qabs, Qsca, Qbk, Qpr, g, Albedo;
    double ti = 0.0, tf = 90.0;
    int nt = 0;
    
    if (argc < 5) {
      printf("Insufficient parameters.\n");
      printf("Usage: %s -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] [-t ti tf nt] [-c comment]\n", argv[0]);
      return -1;
    }
    
    strcpy(comment, "");
    for (i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-l") == 0) {
        i++;
        L = atoi(argv[i]);
        x.resize(L);
        m.resize(L);
        if (argc < 3*(L + 1)) {
          printf("Insufficient parameters.\nUsage: %s -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] [-t ti tf nt] [-c comment]\n", argv[0]);
          return -1;
        } else {
          for (l = 0; l < L; l++) {
            i++;
            x[l] = atof(argv[i]);
            i++;
            m[l] = std::complex<double>(atof(argv[i]), atof(argv[i + 1]));
            i++;
          }
        }
      } else if (strcmp(argv[i], "-t") == 0) {
        i++;
        ti = atof(argv[i]);
        i++;
        tf = atof(argv[i]);
        i++;
        nt = atoi(argv[i]);

        Theta.resize(nt);
        S1.resize(nt);
        S2.resize(nt);
      } else if (strcmp(argv[i], "-c") == 0) {
        i++;
        strcpy(comment, argv[i]);
        has_comment = 1;
      } else { i++; }
    }
    
    if (nt < 0) {
      printf("Error reading Theta.\n");
      return -1;
    } else if (nt == 1) {
      Theta[0] = ti*PI/180.0;
    } else {
      for (i = 0; i < nt; i++) {
      Theta[i] = (ti + (double)i*(tf - ti)/(nt - 1))*PI/180.0;
      }
    }

    nMie(L, x, m, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, S1, S2);

    if (has_comment) {
      printf("%6s, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e\n", comment, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo);
    } else {
      printf("%+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e\n", Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo);
    }
    
    if (nt > 0) {
      printf(" Theta,         S1.r,         S1.i,         S2.r,         S2.i\n");
      
      for (i = 0; i < nt; i++) {
        printf("%6.2f, %+.5e, %+.5e, %+.5e, %+.5e\n", Theta[i]*180.0/PI, S1[i].real(), S1[i].imag(), S2[i].real(), S2[i].imag());
      }
    }
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}


