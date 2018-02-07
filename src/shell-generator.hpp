#ifndef SRC_SHELL_GENERATOR_H_
#define SRC_SHELL_GENERATOR_H_
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
//   @brief  Generates points for integration on sphere surface
#include <complex>
#include <string>
#include <utility>
#include <vector>
namespace shell_generator {
    template<class T> inline T pow2(const T value) {return value*value;}
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  template<typename T, typename A> inline
  T dot(std::vector<T,A> const& a,
        std::vector<T,A> const& b) {
    //return a[0]*std::conj(b[0])+a[1]*std::conj(b[1])+a[2]*std::conj(b[2]);
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    // return std::inner_product(begin(a), end(a), begin(b),
    //                           static_cast<T>(0.0));
  }  
  template<typename T, typename A> inline
  std::vector<T> cross(std::vector<T,A> const& a,
        std::vector<T,A> const& b) {
    std::vector<T> c = {
      a[1]*b[2]-a[2]*b[1],
      a[2]*b[0]-a[0]*b[2],
      a[0]*b[1]-a[1]*b[0]
    };
    return c;
  }  
  // template<typename T, typename A> inline
  // std::vector<std::vector<T> > dyadic(std::vector<T,A> const& a,
  //       std::vector<T,A> const& b) {
  //   std::vector<std::vector<T> > c = {
  //     {a[0]*b[0], a[0]*b[1], a[0]*b[2]},
  //     {a[1]*b[0], a[1]*b[1], a[1]*b[2]},
  //     {a[2]*b[0], a[2]*b[1], a[2]*b[2]}
  //   };
  //   return c;
  // }  
  // template<typename T, typename A> inline
  // std::vector<std::vector<T>,A >
  // operator+(std::vector<std::vector<T>,A > const& a,
  //           std::vector<std::vector<T>,A > const& b) {
  //   std::vector<std::vector<T>,A > c = {
  //     {a[0][0]+b[0][0], a[0][1]+b[0][1], a[0][2]+b[0][2]},
  //     {a[1][0]+b[1][0], a[1][1]+b[1][1], a[1][2]+b[1][2]},
  //     {a[2][0]+b[2][0], a[2][1]+b[2][1], a[2][2]+b[2][2]}
  //   };
  //   return c;
  // }  
  template<typename T, typename A> inline
  std::vector<T,A >
  operator+(std::vector<T,A > const& a,
            std::vector<T,A > const& b) {
    std::vector<T,A > c = {a[0]+b[0],
                           a[1]+b[1],
                           a[2]+b[2]};
    return c;
  }  
  template<typename T, typename A> inline
  std::vector<T,A >
  operator-(std::vector<T,A > const& a,
            std::vector<T,A > const& b) {
    std::vector<T,A > c = {a[0]-b[0],
                           a[1]-b[1],
                           a[2]-b[2]};
    return c;
  }  

  // template<typename T, typename A> inline
  // std::vector<std::vector<T>,A >
  // operator-(std::vector<std::vector<T>,A > const& a,
  //           std::vector<std::vector<T>,A > const& b) {
  //   std::vector<std::vector<T>,A > c = {
  //     {a[0][0]-b[0][0], a[0][1]-b[0][1], a[0][2]-b[0][2]},
  //     {a[1][0]-b[1][0], a[1][1]-b[1][1], a[1][2]-b[1][2]},
  //     {a[2][0]-b[2][0], a[2][1]-b[2][1], a[2][2]-b[2][2]}
  //   };
  //   return c;
  // }  
  // template<typename T, typename A> inline
  // std::vector<std::vector<T>,A >
  // real(std::vector<std::vector<T>,A > const& a) {
  //   std::vector<std::vector<T>,A > c = {
  //     {a[0][0].real(), a[0][1].real(), a[0][2].real()},
  //     {a[1][0].real(), a[1][1].real(), a[1][2].real()},
  //     {a[2][0].real(), a[2][1].real(), a[2][2].real()}
  //   };
  //   return c;
  // }  
  template<typename T, typename A> inline
  std::vector<T>
  real(std::vector<std::complex<T>,A > const& a) {
    std::vector<T> c = {a[0].real(),
                        a[1].real(),
                        a[2].real()};
    return c;
  }
  template<typename T, typename A> inline T
  norm(std::vector<T,A> a){
    T norm_value = 0;
    for (auto coord:a)
      norm_value += pow2(coord);
    return std::sqrt(norm_value);
   }

  template<typename T, typename A> inline T
  real(std::complex<T> const& a) {    
    return a.real();
  }  
  template<typename T, typename A> 
  std::vector< std::complex<T>,A > vconj(std::vector< std::complex<T>,A > const& a) {
    std::vector< std::complex<T>,A > b = {std::conj(a[0]),
                                          std::conj(a[1]),
                                          std::conj(a[2]) };
    // for (auto elem : a)
    //   b.push_back(std::conj(elem));
    return b;
  }  
  template<typename T1, typename T2, typename A> 
  std::vector<T2,A> operator*(T1 const& a, std::vector< T2,A > const& b) {
    std::vector<T2,A > c = {a*b[0],
                            a*b[1],
                            a*b[2]};
    // for (auto elem : b)
    //   c.push_back(a*elem);
    return c;
  }  
  template<typename T1, typename T2, typename A> 
  std::vector<T2,A> operator/(std::vector< T2,A > const& b, T1 const& a) {
    std::vector<T2,A > c = {b[0]/a,
                            b[1]/a,
                            b[2]/a};
    // for (auto elem : b)
    //   c.push_back(a*elem);
    return c;
  }  

  

  class ShellGenerator {  // will throw for any error
   public:
    void SetField(std::vector<std::vector< std::complex<double> > > &E,
                  std::vector<std::vector< std::complex<double> > > &H) {E_ = E; H_=H;};
    void SetFieldSph(std::vector<std::vector< std::complex<double> > > &E,
                  std::vector<std::vector< std::complex<double> > > &H) {Es_ = E; Hs_=H;};
    std::vector<double> IntegrateByFacesQuadrature2();
    std::vector<double> (*ValueAtPoint)(std::vector<double> point) = nullptr;
    std::vector<double> GetChargeField(std::vector<double> point);
    std::vector<double> GetDipoleField(std::vector<double> point);
    std::vector<double> EvaluateDiffForce(const std::vector< std::complex<double> > &E,
                                          const std::vector< std::complex<double> > &H,
                                          const std::vector< std::complex<double> > &sph_unit); 
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
    void PrintVerts();
    void Refine();
    void Rescale(double scale);
    void RotateX(double angle);
    void RotateY(double angle);
    void RotateZ(double angle);
    void MoveX(double delta_x);
    void MoveY(double delta_y);
    void MoveZ(double delta_z);
    void SetInitialVerticesIcosahedron();
    void SetInitialVerticesTetrahedron();

  private:
    double scale_ = 1.0; // Mesh scale (radius of bound sphere)
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
