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
//   This program returns expansion coefficents of Mie series
#include <complex>
#include <cstdio>
#include <string>
#include "../src/nmie-applied.hpp"
#include "../src/nmie-applied-impl.hpp"
#include "../src/nmie-precision.hpp"
#include "./read-spectra.h"

// template<class T> inline T pow2(const T value) {return value*value;}
int main(int argc, char *argv[]) {
  using namespace nmie ;
  try {
    read_spectra::ReadSpectra Si_index, Ag_index;
    read_spectra::ReadSpectra plot_core_index_, plot_TiN_;
    std::string core_filename("Si-int.txt");
    //std::string core_filename("Ag.txt");
    //std::string TiN_filename("TiN.txt");
    std::string TiN_filename("Ag-int.txt");
    //std::string TiN_filename("Si.txt");
    std::string shell_filename(core_filename);

    nmie::MultiLayerMieApplied<nmie::FloatType> multi_layer_mie;  
    const std::complex<double> epsilon_Si(18.4631066585, 0.6259727805);
    const std::complex<double> epsilon_Ag(-8.5014154589, 0.7585845411);
    const std::complex<double> index_Si = std::sqrt(epsilon_Si);
    const std::complex<double> index_Ag = std::sqrt(epsilon_Ag);
    double WL=500; //nm
    double core_width = 5.27; //nm Si
    double inner_width = 8.22; //nm Ag
    double outer_width = 67.91; //nm  Si
    bool isSiAgSi = true;
    double delta_width = 25.0;
    //bool isSiAgSi = false;
    if (isSiAgSi) {
      core_width = 5.27; //nm Si
      inner_width = 8.22; //nm Ag
      outer_width = 67.91; //nm  Si
      multi_layer_mie.AddTargetLayer(core_width, index_Si);
      multi_layer_mie.AddTargetLayer(inner_width, index_Ag);
      multi_layer_mie.AddTargetLayer(outer_width+delta_width, index_Si);
    } else {
      inner_width = 31.93; //nm Ag
      outer_width = 4.06; //nm  Si
      multi_layer_mie.AddTargetLayer(inner_width, index_Ag);
      multi_layer_mie.AddTargetLayer(outer_width+delta_width, index_Si);
    }

    for (int dd = 0; dd<50; ++dd) {
      delta_width = dd;
    FILE *fp;
    std::string fname = "absorb-layered-spectra-d"+std::to_string(dd)+".dat";
    fp = fopen(fname.c_str(), "w");

    multi_layer_mie.SetWavelength(WL);
    multi_layer_mie.RunMieCalculation();
    
    double Qabs = static_cast<double>(multi_layer_mie.GetQabs());
    printf("Qabs = %g\n", Qabs);
    std::vector< std::vector<std::complex<nmie::FloatType> > > aln, bln, cln, dln;
    multi_layer_mie.GetExpanCoeffs(aln, bln, cln, dln);
    std::vector< std::vector<std::complex<double> > > d_aln =
      nmie::ConvertComplexVectorVector<double>(aln);
    std::string str = std::string("#WL ");
    for (int l = 0; l<d_aln.size(); ++l) {
      for (int n = 0; n<3; ++n) {
	str+="|a|^2+|d|^2_ln"+std::to_string(l)+std::to_string(n)+" "
	  + "|b|^2+|c|^2_ln"+std::to_string(l)+std::to_string(n)+" ";
      }
    }
    str+="\n";
    fprintf(fp, "%s", str.c_str());
    fprintf(fp, "# |a|+|d|");
    str.clear();

    double from_WL = 400;
    double to_WL = 600;
    int total_points = 401;
    double delta_WL = std::abs(to_WL - from_WL)/(total_points-1.0);
    Si_index.ReadFromFile(core_filename).ResizeToComplex(from_WL, to_WL, total_points)
      .ToIndex();
    Ag_index.ReadFromFile(TiN_filename).ResizeToComplex(from_WL, to_WL, total_points)
      .ToIndex();
    auto Si_data = Si_index.GetIndex();
    auto Ag_data = Ag_index.GetIndex();
    for (int i=0; i < Si_data.size(); ++i) {
      const double& WL = Si_data[i].first;
      const std::complex<double>& Si = Si_data[i].second;
      const std::complex<double>& Ag = Ag_data[i].second;
      str+=std::to_string(WL);
      multi_layer_mie.ClearTarget();      
      if (isSiAgSi) {
	multi_layer_mie.AddTargetLayer(core_width, Si);
	multi_layer_mie.AddTargetLayer(inner_width, Ag);
	multi_layer_mie.AddTargetLayer(outer_width+delta_width, Si);
      } else {
	inner_width = 31.93; //nm Ag
	outer_width = 4.06; //nm  Si
	multi_layer_mie.AddTargetLayer(inner_width, Ag);
	multi_layer_mie.AddTargetLayer(outer_width+delta_width, Si);
      }
      multi_layer_mie.SetWavelength(WL);
      multi_layer_mie.RunMieCalculation();
      multi_layer_mie.GetQabs();
      multi_layer_mie.GetExpanCoeffs(aln, bln, cln, dln);
      for (int l = 0; l<aln.size(); ++l) {
	for (int n = 0; n<3; ++n) {
	  str+=" "+std::to_string(static_cast<double>(pow2(std::abs(aln[l][n]))+
                                                      pow2(std::abs(dln[l][n]))))
	    + " "
	    + std::to_string(static_cast<double>(pow2(std::abs(bln[l][n]))
                                                 + pow2(std::abs(cln[l][n])) ));

	  // str+=" "+std::to_string(aln[l][n].real() - pow2(std::abs(aln[l][n]))
	  // 			  +dln[l][n].real() - pow2(std::abs(dln[l][n])))
	  //   + " "
	  //   + std::to_string(bln[l][n].real() - pow2(std::abs(bln[l][n]))
	  // 			  +cln[l][n].real() - pow2(std::abs(cln[l][n])) );
	}
      }
      str+="\n";
      fprintf(fp, "%s", str.c_str());
      str.clear();
    }

    fclose(fp);
    }

    // WL = 500;
    // multi_layer_mie.SetWavelength(WL);
    // multi_layer_mie.RunMieCalculation();
    // multi_layer_mie.GetQabs();
    // multi_layer_mie.GetExpanCoeffs(aln, bln, cln, dln);

    // printf("\n Scattering");
    // for (int l = 0; l<aln.size(); ++l) {
    //   int n = 0;
    //   printf("aln[%i][%i] = %g, %gi\n", l, n+1, aln[l][n].real(), aln[l][n].imag());
    //   printf("bln[%i][%i] = %g, %gi\n", l, n+1, bln[l][n].real(), bln[l][n].imag());
    //   printf("cln[%i][%i] = %g, %gi\n", l, n+1, cln[l][n].real(), cln[l][n].imag());
    //   printf("dln[%i][%i] = %g, %gi\n", l, n+1, dln[l][n].real(), dln[l][n].imag());
    //   n = 1;
    //   printf("aln[%i][%i] = %g, %gi\n", l, n+1, aln[l][n].real(), aln[l][n].imag());
    //   printf("bln[%i][%i] = %g, %gi\n", l, n+1, bln[l][n].real(), bln[l][n].imag());
    //   printf("cln[%i][%i] = %g, %gi\n", l, n+1, cln[l][n].real(), cln[l][n].imag());
    //   printf("dln[%i][%i] = %g, %gi\n", l, n+1, dln[l][n].real(), dln[l][n].imag());
    //   // n = 2;
    //   // printf("aln[%i][%i] = %g, %gi\n", l, n+1, aln[l][n].real(), aln[l][n].imag());
    //   // printf("bln[%i][%i] = %g, %gi\n", l, n+1, bln[l][n].real(), bln[l][n].imag());
    //   // printf("cln[%i][%i] = %g, %gi\n", l, n+1, cln[l][n].real(), cln[l][n].imag());
    //   // printf("dln[%i][%i] = %g, %gi\n", l, n+1, dln[l][n].real(), dln[l][n].imag());
    // }
    // printf("\n Absorbtion\n");
    // for (int l = 0; l<aln.size(); ++l) {
    //   if (l == aln.size()-1) printf(" Total ");
    //   printf("===== l=%i   =====\n", l);
    //   int n = 0;
    //   printf("aln[%i][%i] = %g\n", l, n+1, aln[l][n].real() - pow2(std::abs(aln[l][n])));
    //   printf("bln[%i][%i] = %g\n", l, n+1, bln[l][n].real() - pow2(std::abs(bln[l][n])));
    //   printf("cln[%i][%i] = %g\n", l, n+1, cln[l][n].real() - pow2(std::abs(cln[l][n])));
    //   printf("dln[%i][%i] = %g\n", l, n+1, dln[l][n].real() - pow2(std::abs(dln[l][n])));
    //   n = 1;
    //   printf("aln[%i][%i] = %g\n", l, n+1, aln[l][n].real() - pow2(std::abs(aln[l][n])));
    //   printf("bln[%i][%i] = %g\n", l, n+1, bln[l][n].real() - pow2(std::abs(bln[l][n])));
    //   printf("cln[%i][%i] = %g\n", l, n+1, cln[l][n].real() - pow2(std::abs(cln[l][n])));
    //   printf("dln[%i][%i] = %g\n", l, n+1, dln[l][n].real() - pow2(std::abs(dln[l][n])));
    //   // n = 2;
    //   // printf("aln[%i][%i] = %g\n", l, n+1, aln[l][n].real() - pow2(std::abs(aln[l][n])));
    //   // printf("bln[%i][%i] = %g\n", l, n+1, bln[l][n].real() - pow2(std::abs(bln[l][n])));
    //   // printf("cln[%i][%i] = %g\n", l, n+1, cln[l][n].real() - pow2(std::abs(cln[l][n])));
    //   // printf("dln[%i][%i] = %g\n", l, n+1, dln[l][n].real() - pow2(std::abs(dln[l][n])));
    // }


  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}


