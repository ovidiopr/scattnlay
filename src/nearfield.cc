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
#include "nmie.hpp"

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
    std::vector<std::string> args;
    args.assign(argv, argv + argc);
    std::string error_msg(std::string("Insufficient parameters.\nUsage: ") + args[0]
                          + " -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] "
                          + " -p xi xf nx yi yf ny zi zf nz [-c comment]\n");
    enum mode_states {read_L, read_x, read_mr, read_mi, read_xi, read_xf, read_nx, read_yi, read_yf, read_ny, read_zi, read_zf, read_nz, read_comment};
    // for (auto arg : args) std::cout<< arg <<std::endl;
    std::string comment;
    int has_comment = 0;
    unsigned int L = 0;
    std::vector<double> x, Xp, Yp, Zp;
    std::vector<std::complex<double> > m;
    std::vector<std::vector<std::complex<double> > > E, H;

    double xi = 0.0, xf = 0.0, yi = 0.0, yf = 0.0, zi = 0.0, zf = 0.0;
    double dx = 0.0, dy = 0.0, dz = 0.0;
    int nx = 0, ny = 0, nz = 0;
    long total_points = 0;
    if (argc < 5) throw std::invalid_argument(error_msg);

    int mode = -1;
    double tmp_mr;
    for (auto arg : args) {
      // For each arg in args list we detect the change of the current
      // read mode or read the arg. The reading args algorithm works
      // as a finite-state machine.

      // Detecting new read mode (if it is a valid -key)
      if (arg == "-l") {
        mode = read_L;
        continue;
      }

      if (arg == "-p") {
        if ((mode != read_x) && (mode != read_comment))
          throw std::invalid_argument(std::string("Unfinished layer!\n") + error_msg);
        mode = read_xi;
        continue;
      }

      if (arg == "-c") {
        if ((mode != read_x) && (mode != read_nz))
          throw std::invalid_argument(std::string("Unfinished layer or theta!\n") + error_msg);
        mode = read_comment;
        continue;
      }

      // Reading data. For invalid date the exception will be thrown
      // with the std:: and catched in the end.
      if (mode == read_L) {
        L = std::stoi(arg);
        mode = read_x;
        continue;
      }

      if (mode == read_x) {
        x.push_back(std::stod(arg));
        mode = read_mr;
        continue;
      }

      if (mode == read_mr) {
        tmp_mr = std::stod(arg);
        mode = read_mi;
        continue;
      }

      if (mode == read_mi) {
        m.push_back(std::complex<double>( tmp_mr,std::stod(arg) ));
        mode = read_x;
        continue;
      }

      if (mode == read_xi) {
        xi = std::stod(arg);
        mode = read_xf;
        continue;
      }

      if (mode == read_xf) {
        xf = std::stod(arg);
        mode = read_nx;
        continue;
      }

      if (mode == read_nx) {
        nx = std::stoi(arg);
        mode = read_yi;
        continue;
      }

      if (mode == read_yi) {
        yi = std::stod(arg);
        mode = read_yf;
        continue;
      }

      if (mode == read_yf) {
        yf = std::stod(arg);
        mode = read_ny;
        continue;
      }

      if (mode == read_ny) {
        ny = std::stoi(arg);
        mode = read_zi;
        continue;
      }

      if (mode == read_zi) {
        zi = std::stod(arg);
        mode = read_zf;
        continue;
      }

      if (mode == read_zf) {
        zf = std::stod(arg);
        mode = read_nz;
        continue;
      }

      if (mode == read_nz) {
        nz = std::stoi(arg);
        total_points = nx*ny*nz;
        if (total_points <= 0)
          throw std::invalid_argument(std::string("Nothing to do! You must define the grid to calculate the fields.\n") + error_msg);

        Xp.resize(total_points);
        Yp.resize(total_points);
        Zp.resize(total_points);

        E.resize(total_points);
        H.resize(total_points);
        for (long i = 0; i < total_points; i++) {
          E[i].resize(3);
          H[i].resize(3);
        }
        continue;
      }

      if (mode ==  read_comment) {
        comment = arg;
        has_comment = 1;
        continue;
      }
    }

    if ( (x.size() != m.size()) || (L != x.size()) )
      throw std::invalid_argument(std::string("Broken structure!\n") + error_msg);
    if ( (0 == m.size()) || ( 0 == x.size()) )
      throw std::invalid_argument(std::string("Empty structure!\n") + error_msg);

    if (nx == 1)
      dx = 0.0;
    else
      dx = (xf - xi)/(nx - 1);

    if (ny == 1)
      dy = 0.0;
    else
      dy = (yf - yi)/(ny - 1);

    if (nz == 1)
      dz = 0.0;
    else
      dz = (zf - zi)/(nz - 1);

    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          Xp[i*ny*nz + j*nz + k] = xi + (double)i*dx;
          Yp[i*ny*nz + j*nz + k] = yi + (double)j*dy;
          Zp[i*ny*nz + j*nz + k] = zi + (double)k*dz;
        }
      }
    }

    nmie::nField(L, -1, x, m, -1, total_points, Xp, Yp, Zp, E, H);

    if (has_comment)
      printf("%6s\n", comment.c_str());

    if (total_points > 0) {
      printf("         X,          Y,          Z,         Ex.r,         Ex.i,         Ey.r,         Ey.i,         Ez.r,         Ez.i,         Hx.r,         Hx.i,         Hy.r,         Hy.i,         Hz.r,         Hz.i\n");

      for (long i = 0; i < total_points; i++) {
        printf("%10.7f, %10.7f, %10.7f, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e\n",
               Xp[i], Yp[i], Zp[i],
               E[i][0].real(), E[i][0].imag(), E[i][1].real(), E[i][1].imag(), E[i][2].real(), E[i][2].imag(),
               H[i][0].real(), H[i][0].imag(), H[i][1].real(), H[i][1].imag(), H[i][2].real(), H[i][2].imag());
      }
    }


  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }
    return 0;
}


