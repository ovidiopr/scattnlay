//**********************************************************************************//
//    Copyright (C) 2009-2016  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2016  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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
//   This program evaluates forces acting on the nanoparticle under irradiaton.
#include <complex>
#include <cstdio>
#include <string>
#include <iostream>
#include "../src/nmie.hpp"
#include "../src/nmie-impl.hpp"
#include "../src/nmie-applied.hpp"
#include "../src/nmie-applied-impl.hpp"
#include "../src/shell-generator.hpp"
int main(int argc, char *argv[]) {
  try {
    const double pi = 3.1415926535897932384626433832795;
    nmie::MultiLayerMieApplied<double> multi_layer_mie;  
    // const std::complex<double> epsilon_Si(18.4631066585, 0.6259727805);
    //    const std::complex<double> epsilon_Ag(-8.5014154589, 0.7585845411);
    // const std::complex<double> index_Si = std::sqrt(epsilon_Si);
    const std::complex<double> index_Si(2.0);
    //    const std::complex<double> index_Ag = std::sqrt(epsilon_Ag);
    double WL=500; //nm
    // double core_width = 5.27; //nm Si
    // double inner_width = 8.22; //nm Ag
    double outer_width = 67.91; //nm  Si
    // core_width = 5.27; //nm Si
    // inner_width = 8.22; //nm Ag
    //outer_width = WL/(2.0*pi); //nm  Si
    outer_width = 50; //nm  Si
    // multi_layer_mie.AddTargetLayer(core_width, index_Si);
    // multi_layer_mie.AddTargetLayer(inner_width, index_Ag);
    multi_layer_mie.AddTargetLayer(outer_width, index_Si);
    multi_layer_mie.SetWavelength(WL);
    multi_layer_mie.RunMieCalculation();
    double Qsca = multi_layer_mie.GetQsca();
    printf("Qsca = %g\n", Qsca);

    shell_generator::ShellGenerator shell;
    shell.Init();
    // shell.Refine();
    // shell.Refine();
    //shell.Refine();

    //double scale = 2.0*pi*(core_width + inner_width +  outer_width)/WL*8.0001;  //Integration sphere radius.
    double scale = 2.0*pi*(outer_width)/WL*1.0001;  //Integration sphere radius.
    //double scale = 1.0001;  //Integration sphere radius.
    
    shell.Rescale(scale);
    // shell.RotateX(pi/2.0);
    // shell.RotateY(pi/2.0);
    // shell.RotateZ(pi/2.0);
    //shell.PrintVerts();
    auto points = shell.GetVerticesT();
    multi_layer_mie.SetFieldPointsSP(points);
    multi_layer_mie.RunFieldCalculation();
    auto E = nmie::ConvertComplexVectorVector<double>(multi_layer_mie.GetFieldE());
    auto H = nmie::ConvertComplexVectorVector<double>(multi_layer_mie.GetFieldH());
    shell.SetField(E,H);
    auto F = shell.Integrate();
    std::cout<<"F: " <<F[0]<<", "<< F[1] <<", "<<F[2] << std::endl<< std::endl;
    F = shell.IntegrateByComp();
    std::cout<<"F: " <<F[0]<<", "<< F[1] <<", "<<F[2] << std::endl;
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}


