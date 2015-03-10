/**
 * @file   test-negative-epsilon.cc
 * @author Konstantin Ladutenko <kostyfisik at gmail (.) com>
 * @date   Mon Mar  9 13:21:37 2015
 * 
 * @brief  test negative epsilon case
 * 
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
 * 
 */
#include "nmie.h"
#include <stdio.h>
template<class T> inline T pow2(const T value) {return value*value;};
const double PI=3.14159265358979323846;  
nmie::MultiLayerMie multi_layer_mie_;  
double lambda_work_ = 3.75; // cm
double a_ = 0.75*lambda_work_;  // 2.8125 cm - size of PEC core
double min_index_ = 1e-11;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateScatterOnlyIndex(std::vector<double> input) {
  double Qsca;
  std::vector<std::complex<double>> cindex;
  cindex.clear();
  double k = min_index_, n=min_index_;
  for (double epsilon : input) {
    // sqrt(epsilon) = n + i*k
    k = min_index_;
    n=min_index_;
    if (epsilon > 0.0) n=std::sqrt(epsilon);
    else k = std::sqrt(-epsilon);
    if (n < min_index_) n = min_index_;
    if (k < min_index_) k = min_index_;
    //printf("eps= %g, n=%g, k=%g\n", epsilon, n, k);
    cindex.push_back(std::complex<double>(n, k));
  }
  multi_layer_mie_.SetCoatingIndex(cindex);
  try {
    multi_layer_mie_.RunMieCalculations();
    Qsca = multi_layer_mie_.GetQsca();
  } catch( const std::invalid_argument& ia ) {
    Qsca = 0;
    printf("#");
    // Will catch if  multi_layer_mie_ fails or other errors.
    //std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  double total_r = multi_layer_mie_.GetTotalRadius();
  return Qsca*PI*pow2(total_r);
}

int main(int argc, char *argv[]) {
  try {
    // Only PEC target
    multi_layer_mie_.SetTargetPEC(a_);
    multi_layer_mie_.SetWavelength(lambda_work_);
    multi_layer_mie_.RunMieCalculations();
    double PEC_Qsca = multi_layer_mie_.GetQsca();
    double PEC_r = multi_layer_mie_.GetTotalRadius();
    double PEC_RCS = PEC_Qsca*PI*pow2(PEC_r);
    
    // PEC target covered with with air layer
    multi_layer_mie_.SetCoatingWidth({0.1});
    multi_layer_mie_.SetCoatingIndex({{1.0,0.0}});
    multi_layer_mie_.RunMieCalculations();
    double Qsca1 = multi_layer_mie_.GetQsca();
    double total_r1 = multi_layer_mie_.GetTotalRadius();
    double initial_RCS1 = Qsca1*PI*pow2(total_r1);
    printf("RCS = %g cm^2 with (r=%g) and  RCS=%g cm^2 without (r=%g)air coating.\n",
	   initial_RCS1, total_r1, 
	   PEC_RCS, PEC_r);

    //multi_layer_mie.SetMaxTermsNumber(150);

    // Bi-layer, inner layer =  lossless metall.
    std::vector<std::complex<double>> cindex;
    cindex.clear();
    double n=min_index_;
    double k=std::sqrt(0.29);
    cindex.push_back(std::complex<double>(n, k));
    n=std::sqrt(24.6);
    k=min_index_;
    cindex.push_back(std::complex<double>(n, k));

    multi_layer_mie_.SetCoatingWidth({0.1,0.1});
    multi_layer_mie_.SetCoatingIndex(cindex);
    multi_layer_mie_.RunMieCalculations();
    double Qsca = multi_layer_mie_.GetQsca();
    double total_r = multi_layer_mie_.GetTotalRadius();
    double initial_RCS = Qsca*PI*pow2(total_r);
    printf("RCS=%g for bi-layer coating (total R=%g).\n", initial_RCS,total_r);

    n=1.0;
    k=min_index_;
    cindex.push_back(std::complex<double>(n, k));
    multi_layer_mie_.SetCoatingWidth({0.1,0.1,0.1});
    multi_layer_mie_.SetCoatingIndex(cindex);
    multi_layer_mie_.RunMieCalculations();
    Qsca = multi_layer_mie_.GetQsca();
    total_r = multi_layer_mie_.GetTotalRadius();
    initial_RCS = Qsca*PI*pow2(total_r);
    printf("RCS=%g for bi-layer+air coating (total R=%g).\n", initial_RCS,total_r);


    //multi_layer_mie_.SetMaxTermsNumber(15);
    multi_layer_mie_.SetCoatingWidth({0.1,0.1});
    printf("With %g  coating= (26.24).\n",
	   EvaluateScatterOnlyIndex({-0.29, 24.6}));
    multi_layer_mie_.SetCoatingWidth({0.1,0.1,0.1});
    printf("With %g  coating> (26.24).\n",
	   EvaluateScatterOnlyIndex({-0.29, 24.6, 1.0}));
    //multi_layer_mie_.SetMaxTermsNumber(-1);

      // 26.24: 25||   -0.29   +24.60 
      // 28.48: 38||   -0.29   +24.60    +1.00
    
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}
