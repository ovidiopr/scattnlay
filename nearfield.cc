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

using namespace std;

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
    vector<string> args;
    args.assign(argv, argv + argc);
    string error_msg(string("Insufficient parameters.\nUsage: ") + args[0]
                          + " -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] "
                          + "[-p xi xf nx yi yf ny zi zf nz] [-c comment]\n");
    enum mode_states {read_L, read_x, read_mr, read_mi, read_xi, read_xf, read_nx, read_yi, read_yf, read_ny, read_zi, read_zf, read_nz, read_comment};
    // for (auto arg : args) cout<< arg <<endl;
    string comment;
    int has_comment = 0;
    unsigned int L = 0;
    vector<double> x, Xp, Yp, Zp;
    vector<complex<double> > m;
    vector<vector<complex<double> > > E, H;

    double xi = 0.0, xf = 0.0, yi = 0.0, yf = 0.0, zi = 0.0, zf = 0.0;
    double dx = 0.0, dy = 0.0, dz = 0.0;
    int nx = 0, ny = 0, nz = 0;
    long total_points = 0;
    if (argc < 5) throw invalid_argument(error_msg);

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
          throw invalid_argument(string("Unfinished layer!\n") + error_msg);
        mode = read_xi;
        continue;
      }

      if (arg == "-c") {
        if ((mode != read_x) && (mode != read_nz))
          throw invalid_argument(string("Unfinished layer or theta!\n") + error_msg);
        mode = read_comment;
        continue;
      }

      // Reading data. For invalid date the exception will be thrown
      // with the  and catched in the end.
      if (mode == read_L) {
        L = stoi(arg);
        mode = read_x;
        continue;
      }

      if (mode == read_x) {
        x.push_back(stod(arg));
        mode = read_mr;
        continue;
      }

      if (mode == read_mr) {
        tmp_mr = stod(arg);
        mode = read_mi;
        continue;
      }

      if (mode == read_mi) {
        m.push_back(complex<double>( tmp_mr,stod(arg) ));
        mode = read_x;
        continue;
      }

      if (mode == read_xi) {
        xi = stod(arg);
        mode = read_xf;
        continue;
      }

      if (mode == read_xf) {
        xf = stod(arg);
        mode = read_nx;
        continue;
      }

      if (mode == read_nx) {
        nx = stoi(arg);
        mode = read_yi;
        continue;
      }

      if (mode == read_yi) {
        yi = stod(arg);
        mode = read_yf;
        continue;
      }

      if (mode == read_yf) {
        yf = stod(arg);
        mode = read_ny;
        continue;
      }

      if (mode == read_ny) {
        ny = stoi(arg);
        mode = read_zi;
        continue;
      }

      if (mode == read_zi) {
        zi = stod(arg);
        mode = read_zf;
        continue;
      }

      if (mode == read_zf) {
        zf = stod(arg);
        mode = read_nz;
        continue;
      }

      if (mode == read_nz) {
        nz = stoi(arg);
        total_points = nx*ny*nz;
        if (total_points <= 0)
          throw invalid_argument(string("Nothing to do! You must define the grid to calculate the fields.\n") + error_msg);

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
      throw invalid_argument(string("Broken structure!\n") + error_msg);
    if ( (0 == m.size()) || ( 0 == x.size()) )
      throw invalid_argument(string("Empty structure!\n") + error_msg);

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
          Xp[i*ny + j*nz + k] = xi + (double)i*dx;
          Yp[i*ny + j*nz + k] = yi + (double)j*dy;
          Zp[i*ny + j*nz + k] = zi + (double)k*dz;
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


  } catch( const invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    cerr << "Invalid argument: " << ia.what() << endl;
    return -1;
  }
    return 0;
}


