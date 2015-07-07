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
#include "../src/nmie-applied.h"
int main(int argc, char *argv[]) {
  try {
    nmie::MultiLayerMieApplied multi_layer_mie;  
    const std::complex<double> epsilon_Si(18.4631066585, 0.6259727805);
    const std::complex<double> epsilon_Ag(-8.5014154589, 0.7585845411);
    const std::complex<double> index_Si = std::sqrt(epsilon_Si);
    const std::complex<double> index_Ag = std::sqrt(epsilon_Ag);
    const double WL=500; //nm
    double core_width = 5.27; //nm Si
    double inner_width = 8.22; //nm Ag
    double outer_width = 67.91; //nm  Si
    //bool isSiAgSi = true;
    bool isSiAgSi = false;
    if (isSiAgSi) {
      multi_layer_mie.AddTargetLayer(core_width, index_Si);
      multi_layer_mie.AddTargetLayer(inner_width, index_Ag);
      multi_layer_mie.AddTargetLayer(outer_width, index_Si);
    } else {
      inner_width = 31.93; //nm Ag
      outer_width = 4.06; //nm  Si
      multi_layer_mie.AddTargetLayer(inner_width, index_Ag);
      multi_layer_mie.AddTargetLayer(outer_width, index_Si);
    }
    multi_layer_mie.SetWavelength(WL);
    multi_layer_mie.RunMieCalculation();
    double Qabs = multi_layer_mie.GetQabs();
    printf("Qabs = %g\n", Qabs);
    std::vector< std::vector<std::complex<double> > > aln, bln, cln, dln;
    multi_layer_mie.GetExpanCoeffs(aln, bln, cln, dln);
    for (int l = 0; l<aln.size(); ++l) {
    int n = 0;
    printf("aln[%i][%i] = %g, %gi)\n", l, n, aln[l][n].real(), aln[l][n].imag());
    printf("bln[%i][%i] = %g, %gi)\n", l, n, bln[l][n].real(), bln[l][n].imag());
    printf("cln[%i][%i] = %g, %gi)\n", l, n, cln[l][n].real(), cln[l][n].imag());
    printf("dln[%i][%i] = %g, %gi)\n", l, n, dln[l][n].real(), dln[l][n].imag());
  }



  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}


