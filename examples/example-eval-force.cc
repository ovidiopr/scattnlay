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
double scale_ = 1.0;
  const double pi = 3.1415926535897932384626433832795;
//const double pi = 3.1415926535897932384626433832795;
double WL=545; //nm
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
std::vector<double>
EvaluateDiffForce (const std::vector< std::complex<double> > &E,
                   const std::vector< std::complex<double> > &H,
                   const std::vector<std::complex<double> > unit) {
  using namespace shell_generator;
  std::vector<double> P = (1/(2.0))
    *real(
          dot(unit,E)*vconj(E) +
          dot(unit,H)*vconj(H) +
          (-1.0/2.0)*(dot(E,vconj(E))
                      +dot(H,vconj(H))
                      )*unit
          );
  return P;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
std::vector<double> GetChargeField (std::vector<double> point) {
  using namespace shell_generator;
  
  std::vector< std::complex<double> > zero (3,std::complex<double>(0.0,0.0));
  std::vector< std::complex<double> > E = zero;
  std::vector< std::complex<double> > H = zero;
  
  //double charge = 3.14;
  double charge = 1.0;
  double shift = 10;//*(2.0*pi)/WL;
  if (norm(point) < shift) std::cout<<"<";
    std::vector<double> v_shift = {shift, 0.0, 0.0};

    double r = norm(point-v_shift);
    std::vector<std::complex<double> > sph_unit = { point[0]/r, point[1]/r, point[2]/r};
    std::vector<std::complex<double> > unit = { (point[0]-shift)/r, point[1]/r, point[2]/r};

    const double pi = 3.1415926535897932384626433832795;
  //const double pi = 3.1415926535897932384626433832795;
    double ampl = charge/(4.0*pi*pow2(r));      
    E = ampl*unit;
    std::cout << "%% " << real(E[0]) << " "
              << real(E[1]) << " "
              << real(E[2]) << " "
              << std::endl;
    //return EvaluateDiffForce(E,H,sph_unit);

    std::vector< double > gauss (3, 0.0);
    std::complex< double > gauss_value = dot(sph_unit,E);
    gauss[0] = real(gauss_value);
    return gauss;


    // // Test Poynting vector integration
    // std::vector<double> unit = { vert[0]/r, vert[1]/r, vert[2]/r};
    // std::vector<double> P = (1/(2.0))
    //   *real(cross(E,vconj(H)));
    //integral[0] = integral[0] + per_face_area_[i]*dot(P,unit);

    
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //

int main(int argc, char *argv[]) {
  try {
    nmie::MultiLayerMieApplied<double> multi_layer_mie;  
    const std::complex<double> epsilon_Si(18.4631066585, 0.6259727805);
    //    const std::complex<double> epsilon_Ag(-8.5014154589, 0.7585845411);
    const std::complex<double> index_Si = std::sqrt(epsilon_Si);
    //const std::complex<double> index_Si(3.1,0.00);
    //    const std::complex<double> index_Ag = std::sqrt(epsilon_Ag);
    //double WL=400; //nm
    //double outer_width = 67.91; //nm  Si
    //double outer_width = 40; //nm  Si
    double outer_width = 1; //nm  Si
    //auto shift = 0.0;
    for (int refines=0; refines<1; ++refines) {
      //shell.Refine();
      //for (int i=1; i<165; ++i) {
      for (int i=4; i<5; ++i) {
        //outer_width = 40 + 5*i;
        auto integration_radius = outer_width  + 5*i ;
        //auto integration_radius = 1.0 ;
        //outer_width = 10; //+10*i; //nm  Si
        multi_layer_mie.ClearAllDesign();
        multi_layer_mie.AddTargetLayer(outer_width, index_Si);
        multi_layer_mie.SetWavelength(WL);
        multi_layer_mie.RunMieCalculation();
        //double Qsca = multi_layer_mie.GetQsca();
        //printf("Qsca = %g\t", Qsca);
        scale_ = // 2.0*pi*
          (integration_radius);///WL;//*1.00001;  //Integration sphere radius.
        shell_generator::ShellGenerator shell;
        shell.Init();
        shell.Refine();    //     shell.Refine();         shell.Refine();
        shell.Rescale(scale_);
        auto points = shell.GetVerticesT();
        //auto points = shell.GetFaceCentersT();
        multi_layer_mie.SetFieldPointsSP(points);
        multi_layer_mie.RunFieldCalculation();
        auto E = nmie::ConvertComplexVectorVector<double>(multi_layer_mie.GetFieldE());
        auto H = nmie::ConvertComplexVectorVector<double>(multi_layer_mie.GetFieldH());
        // auto Es = nmie::ConvertComplexVectorVector<double>(multi_layer_mie.GetFieldEs());
        // auto Hs = nmie::ConvertComplexVectorVector<double>(multi_layer_mie.GetFieldHs());
        shell.SetField(E,H);
        //shell.SetFieldSph(Es,Hs);
        //auto F = shell.Integrate();

        shell.ValueAtPoint = &GetChargeField;

        auto F = shell.IntegrateByFacesQuadrature2();
        //auto F = shell.IntegrateByComp();

        std::cout //<< "integrate_R:\t"
          <<std::setprecision(16)
          << (scale_)//*WL/(2.0*pi)
          // << " $$ "<< shell_generator::norm(points[0])
          ;
        
        std::cout<<"\t"
                 <<F[0]<<"\t"<< F[1] <<"\t" <<F[2] << std::endl;

        // auto F1 = shell.IntegrateByComp();
        // std::cout<<"F: " <<F1[0]<<", "<< F1[1] <<", "<<F1[2] << std::endl;        

        // auto F = shell.IntegrateGaussSimple(3.14, 2*pi*outer_width/WL);
        // std::cout<<"\tr: "<<outer_width/2.0<<"\tF: " <<F<< std::endl;        
        
      }

    }  // end for refines
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}


