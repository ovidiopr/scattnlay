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
//sudo aptitude install libgoogle-perftools-dev
//#include <google/heap-profiler.h>
#include "../../src/nmie-applied.h"

timespec diff(timespec start, timespec end);
const double PI=3.14159265358979323846;
template<class T> inline T pow2(const T value) {return value*value;}


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
			  + "[-t ti tf nt] [-c comment]\n");
    enum mode_states {read_L, read_x, read_mr, read_mi, read_ti, read_tf, read_nt, read_comment};
    // for (auto arg : args) std::cout<< arg <<std::endl;
    std::string comment;
    int has_comment = 0;
    int i, l, L = 0;
    std::vector<double> x, Theta;
    std::vector<std::complex<double> > m, S1, S2;
    double Qext, Qabs, Qsca, Qbk, Qpr, g, Albedo;

    std::vector<std::complex<double> > mw, S1w, S2w;
    double Qextw, Qabsw, Qscaw, Qbkw, Qprw, gw, Albedow;

    double ti = 0.0, tf = 90.0;
    int nt = 0;    
    if (argc < 5) throw std::invalid_argument(error_msg);
    
    //strcpy(comment, "");
    // for (i = 1; i < argc; i++) {
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
      if (arg == "-t") {
	if ((mode != read_x) && (mode != read_comment))
	  throw std::invalid_argument(std::string("Unfinished layer!\n")
							 +error_msg);
	mode = read_ti;
	continue;
      }
      if (arg == "-c") {
	if ((mode != read_x) && (mode != read_nt))
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
      if (mode == read_ti) {
	ti = std::stod(arg);
	mode = read_tf;
	continue;
      }
      if (mode == read_tf) {
	tf = std::stod(arg);
	mode = read_nt;
	continue;
      }
      if (mode == read_nt) {
	nt = std::stoi(arg);
        Theta.resize(nt);
        S1.resize(nt);
        S2.resize(nt);
        S1w.resize(nt);
        S2w.resize(nt);
	continue;
      }
      if (mode ==  read_comment) {
	comment = arg;
        has_comment = 1;
	continue;
      }
    }
    if ( (x.size() != m.size()) || (L != x.size()) ) 
      throw std::invalid_argument(std::string("Broken structure!\n")
							 +error_msg);
    if ( (0 == m.size()) || ( 0 == x.size()) ) 
      throw std::invalid_argument(std::string("Empty structure!\n")
							 +error_msg);
    
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
    timespec time1, time2;
    long cpptime_nsec, best_cpp;
    long ctime_nsec, best_c;
    long cpptime_sec, ctime_sec;
    long repeats = 150;
    //HeapProfilerStart("heapprof");    
    do {
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
      for (int i = 0; i<repeats; ++i) {
    	nmie::nMieApplied(L, x, m, nt, Theta, &Qextw, &Qscaw,
    			   &Qabsw, &Qbkw, &Qprw, &gw, &Albedow, S1w, S2w);
      }
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
      cpptime_nsec = diff(time1,time2).tv_nsec;
      cpptime_sec = diff(time1,time2).tv_sec;

      printf("-- C++ time consumed %lg sec\n", (cpptime_nsec/1e9));
      repeats *= 10;
    } while (cpptime_nsec < 1e8 && ctime_nsec < 1e8);

        printf("\n");
    
    if (has_comment) {
      printf("%6s, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e  \n", comment.c_str(), Qextw, Qscaw, Qabsw, Qbkw, Qprw, gw, Albedow);
    } else {
      printf("%+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e  \n", Qextw, Qscaw, Qabsw, Qbkw, Qprw, gw, Albedow);
    }
    
    if (nt > 0) {
      printf(" Theta,         S1.r,         S1.i,         S2.r,         S2.i\n");
      
      for (i = 0; i < nt; i++) {
        printf("%6.2f, %+.5e, %+.5e, %+.5e, %+.5e  \n", Theta[i]*180.0/PI, S1w[i].real(), S1w[i].imag(), S2w[i].real(), S2w[i].imag());
      }
    }

  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}



timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}
