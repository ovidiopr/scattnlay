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
#include "../src/shell-generator.hpp"
int main(int argc, char *argv[]) {
  try {
    shell_generator::ShellGenerator shell;
    shell.Init();
    shell.Refine();
    //shell.Refine();
    //shell.Refine();
    double scale = 1.9;  //Integration sphere radius.
    double shift = 1.1;  //Shift of the charge relative to the sphere center.
    double charge = 11.0; //Coulomb charge.
    shell.Rescale(scale);
    shell.PrintVerts();
    auto points = shell.GetVertices();
    double charge_s = shell.IntegrateGaussSimple(charge, shift);
    std::cout << "Accuracy ( 1==ideal ): " << charge_s/charge << std::endl; 
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}


