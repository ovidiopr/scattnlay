#ifndef SRC_SHELL_GENERATOR_H_
#define SRC_SHELL_GENERATOR_H_
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
//   @brief  Generates points for integration on sphere surface
#include <complex>
#include <string>
#include <utility>
#include <vector>
namespace shell_generator {
  class ShellGenerator {  // will throw for any error
   public:
    /* ShellGenerator& ReadFromFile(std::string filename); */
    /* ShellGenerator& ResizeToComplex(double from_wl, double to_wl, int samples); */
    /* ShellGenerator& ToIndex(); */
    double dist(std::vector<double> a, std::vector<double> b);
    double norm(std::vector<double> a);
    std::vector< std::vector<double> > GetVertices(){return vertices_;};
    std::vector< std::vector<double> > GetVerticesT();
    std::vector< std::vector<double> > GetFaceCentersT();
    std::vector< std::vector<long unsigned int> > GetEdges(){return edges_;};
    std::vector< std::vector<long unsigned int> > GetFaces(){return faces_;};
    void EvalEdgeLength();
    void EvalFaces();
    void GenerateEdgesInit();
    void GenerateFacesInit();
    void Init();
    std::vector<double> Integrate();
    std::vector<double> IntegrateByFaces();
    std::vector<double> IntegrateByComp();
    double IntegrateGauss(double charge, double dist);
    double IntegrateGaussSimple(double charge, double dist);
    void PrintVerts();
    void Refine();
    void Rescale(double scale);
    void RotateX(double angle);
    void RotateY(double angle);
    void RotateZ(double angle);
    void MoveX(double delta_x);
    void MoveY(double delta_y);
    void MoveZ(double delta_z);

    void SetField(std::vector<std::vector< std::complex<double> > > &E,
                  std::vector<std::vector< std::complex<double> > > &H) {E_ = E; H_=H;};
    void SetFieldSph(std::vector<std::vector< std::complex<double> > > &E,
                  std::vector<std::vector< std::complex<double> > > &H) {Es_ = E; Hs_=H;};
    void SetInitialVerticesIcosahedron();
    void SetInitialVerticesTetrahedron();
  private:
    std::vector<std::vector< std::complex<double> > > E_, H_, Es_, Hs_;
    double edge_length_ = 0.0;
    double face_area_ = 0.0;
    std::vector<double> per_face_area_;
    double per_vertice_area_ = 0.0;
    std::vector< std::vector<double> > vertices_, face_centers_;
    std::vector< double > face_surface_; 
    std::vector< std::vector<long unsigned int> > edges_, refined_edges_;
    std::vector< std::vector<long unsigned int> > faces_, refined_faces_;
    // std::vector< std::pair< double, std::complex<double> > > data_complex_;
    // std::vector< std::pair< double, std::complex<double> > > data_complex_index_;
    // void PermittivityToIndex();
  };  // end of class ShellGenerator
}  // end of namespase read_spectra
#endif  // SRC_SHELL_GENERATOR_H_
