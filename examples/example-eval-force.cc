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
    const std::complex<double> index_Si(1.1,0.0);
    //    const std::complex<double> index_Ag = std::sqrt(epsilon_Ag);
    double WL=545; //nm
    //double WL=400; //nm
    //double outer_width = 67.91; //nm  Si
    double outer_width = 4*2*2; //nm  Si
    auto shift = 0.0;
    shell_generator::ShellGenerator shell;
    shell.Init();
    shell.Refine();
    shell.Refine();
    for (int refines=0; refines<1; ++refines) {
      shell.Refine();
      for (int i=0; i<170; ++i) {
        //outer_width = 40 + 5*i;
        auto integration_radius = outer_width  + 5*i ;
        //outer_width = 10; //+10*i; //nm  Si
        multi_layer_mie.ClearAllDesign();
        multi_layer_mie.AddTargetLayer(outer_width, index_Si);
        multi_layer_mie.SetWavelength(WL);
        multi_layer_mie.RunMieCalculation();
        double Qsca = multi_layer_mie.GetQsca();
        //printf("Qsca = %g\t", Qsca);
        double scale = 2.0*pi*(integration_radius)/WL*1.00001;  //Integration sphere radius.
        shell.Rescale(scale);
        // auto points = shell.GetVerticesT();
        auto points = shell.GetFaceCentersT();
        multi_layer_mie.SetFieldPointsSP(points);
        multi_layer_mie.RunFieldCalculation();
        auto E = nmie::ConvertComplexVectorVector<double>(multi_layer_mie.GetFieldE());
        auto H = nmie::ConvertComplexVectorVector<double>(multi_layer_mie.GetFieldH());
        // auto Es = nmie::ConvertComplexVectorVector<double>(multi_layer_mie.GetFieldEs());
        // auto Hs = nmie::ConvertComplexVectorVector<double>(multi_layer_mie.GetFieldHs());
        shell.SetField(E,H);
        //shell.SetFieldSph(Es,Hs);
        // auto F = shell.Integrate();
        //auto F = shell.IntegrateByFaces();
        auto F = shell.IntegrateByComp();
        std::cout << "integrate_R:\t" << scale*WL/(2.0*pi);
        std::cout<<"\tforce:\t" <<F[0]<<"\t"<< F[1] <<"\t"<<F[2] << std::endl;
        // auto F1 = shell.IntegrateByComp();
        // std::cout<<"F: " <<F1[0]<<", "<< F1[1] <<", "<<F1[2] << std::endl;        
      }

    }  // end for refines
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}


