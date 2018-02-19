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

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "shell-generator.hpp"
namespace shell_generator {
  template<class T> inline T pow2(const T value) {return value*value;}
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  template<typename T, typename A> inline
  T dot(std::vector<T,A> const& a,
        std::vector<T,A> const& b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
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
  template<typename T, typename A> inline
  T  real(std::complex<T> const& a) {    
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

  
  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::RotateZ(double angle) {
    for(auto& p : vertices_) {
      double x,y;
      x = std::cos(angle)*p[0]-std::sin(angle)*p[1];
      y = std::sin(angle)*p[0]+std::cos(angle)*p[1];
      p[0] = x;
      p[1] = y;
    }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::RotateY(double angle) {
    for(auto& p : vertices_) {
      double x,z;
      x = std::cos(angle)*p[0]+std::sin(angle)*p[2];
      z = -std::sin(angle)*p[0]+std::cos(angle)*p[2];
      p[0] = x;
      p[2] = z;
    }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::RotateX(double angle) {
    for(auto& p : vertices_) {
      double y,z;
      y = std::cos(angle)*p[1]-std::sin(angle)*p[2];
      z = std::sin(angle)*p[1]+std::cos(angle)*p[2];
      p[1] = y;
      p[2] = z;
    }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::MoveX(double delta) {
    for(auto& p : vertices_) {
      p[0] += delta;
    }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::MoveY(double delta) {
    for(auto& p : vertices_) {
      p[1] += delta;
    }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::MoveZ(double delta) {
    for(auto& p : vertices_) {
      p[2] += delta;
    }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double ShellGenerator::IntegrateGaussSimple(double charge, double shift) {
    double integral = 0.0;
    // return field at point p from the charge, located at (shift, 0,0)
    auto field = [](double charge, double shift, std::vector<double> p){
      double r = std::sqrt(pow2(p[0]-shift) + pow2(p[1]) + pow2(p[2]) );
      //std::cout << "r: " << r << std::endl;
      const double pi = 3.1415926535897932384626433832795;
      double ampl = charge/(4.0*pi*pow2(r));
      std::vector<double> field = {ampl*(p[0]-shift)/r, ampl*(p[1])/r, ampl*(p[2])/r};
      return field;
    };
    //simple 
    // for (auto vert :vertices_) {
    for (long unsigned int i = 0; i<face_centers_.size(); ++i) {
      auto vert = face_centers_[i];
      auto E0 = field(charge, shift, vert);
      // std::cout << "E0: ";
      // for (auto component : E0) std::cout << component << " ";
      // std::cout << std::endl;
      // Vector to unit product
      double r = norm(vert);
      std::vector<double> unit = { vert[0]/r, vert[1]/r, vert[2]/r};
      // std::cout << norm(unit) << std::endl;
      for (int j =0; j < 3; ++j) 
        integral += per_face_area_[i]*unit[j]*E0[j];
    }
    return integral;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double ShellGenerator::IntegrateGauss(double charge, double shift) {
    if (faces_.size() == 0)
      throw std::invalid_argument("Error! Faces were not initialized!\nSee IntegrateGaussSimple for using vertices information only.");
    double integral = 0.0;
    // return field at point p from the charge, located at (shift, 0,0)
    auto field = [](double charge, double shift, std::vector<double> p){
      double r = std::sqrt(pow2(p[0]-shift) + pow2(p[1]) + pow2(p[2]) );
      const double pi = 3.1415926535897932384626433832795;
      double ampl = charge/(4.0*pi*pow2(r));      
      std::vector<double> field = {ampl*(p[0]-shift)/r, ampl*(p[1])/r, ampl*(p[2])/r};
      return field;
    };
    
    for(const auto face : faces_){
      std::vector<double> mean_vector = {0.0, 0.0, 0.0}, mean_point = {0.0, 0.0, 0.0};
      //Get mean
      for (int vert = 0; vert<3; ++vert) { //vertice
        auto point = vertices_[face[vert]];
        auto E0 = field(charge, shift, point);
        for (int i=0; i<3; ++i) {
          mean_vector[i] += E0[i]/3.0;
          mean_point[i] += point[i]/3.0;
        }
      }
      // Vector to unit product
      double r = norm(mean_point);
      std::vector<double> unit = { mean_point[0]/r, mean_point[1]/r, mean_point[2]/r};
      // std::cout << norm(unit) << std::endl;
      for (int i =0; i < 3; ++i) 
        integral += face_area_*
          unit[i]*mean_vector[i];
    }
    return integral;
  }

  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  
  std::vector<double> ShellGenerator::IntegrateByComp() {
    std::vector<double> integral = {0.0, 0.0, 0.0};
    for (unsigned int i=0; i<E_.size(); ++i) {
      auto E = E_[i];
      auto H = H_[i];
      auto vert = face_centers_[i];
      double r = norm(vert);
      std::vector<double> n = { vert[0]/r, vert[1]/r, vert[2]/r};
      const std::vector< std::vector<double> > d = {{1.0, 0.0, 0.0},
                                                    {0.0, 1.0, 0.0},
                                                    {0.0, 0.0, 1.0}};
      
      std::vector<double> F = {0.0, 0.0, 0.0};
      std::complex<double> S(0.0);
      for (int ii = 0; ii < 3; ++ii)
        S += E[ii]*std::conj(E[ii]) + H[ii]*std::conj(H[ii]);
      std::vector< std::vector<std::complex<double> > >
        T = {{0.0, 0.0, 0.0},
             {0.0, 0.0, 0.0},
             {0.0, 0.0, 0.0}};
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          T[i][j] =  E[i]*std::conj(E[j]) + H[i]*std::conj(H[j])
            -1.0/2.0*S*d[i][j];
          F[i] += (1/(2.0/* *4.0*pi */))*real(T[i][j]*n[j]);
        }
      }
      integral = integral + per_face_area_[i]*F;      
    }
    return integral;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  
  std::vector<double> ShellGenerator::IntegrateByFaces() {
    std::vector<double> integral = {0.0, 0.0, 0.0};
    //simple 
    for (long unsigned int i=0; i<E_.size(); ++i) {
      //std::cout << i << " ";
      auto E = E_[i];
      //auto H = 377.0*H_[i];
      auto H = H_[i];
      // auto Es = Es_[i];
      // auto Hs = Hs_[i];
      
      auto vert = face_centers_[i];
      // Vector to unit product
      double r = norm(vert);
      //std::vector<std::complex<double> > unit = { vert[0]/r, vert[1]/r, vert[2]/r};
      // std::cout << norm(unit) << std::endl;
      //const double pi = 3.1415926535897932384626433832795;
      // std::vector<double> P = (1/(2.0))
      //   *real(
      //         dot(unit,E)*vconj(E) +
      //         dot(unit,H)*vconj(H) +
      //         (-1.0/2.0)*(dot(E,vconj(E))
      //                     +dot(H,vconj(H))
      //                     )*unit
      //         );

      // std::vector<double> P = (1/(2.0))
      //   *real(
      //         Es[0]*vconj(E) +
      //         Hs[0]*vconj(H) +
      //         (-1.0/2.0)*(dot(E,vconj(E))
      //                     +dot(H,vconj(H))
      //                     )*unit
      //         );
      // std::vector<double> P = (1/(2.0))
      //   *(
      //     real(dot(unit,E)*vconj(E)) +
      //     real(dot(unit,H)*vconj(H))) +
      //   (-1.0/4.0)*(dot(E,vconj(E))*unit
      //               +dot(H,vconj(H))*unit
      //               )
              
      //     );
      // auto
      // std::cout <<"E "<<E[0]<<", "<< E[1] <<", "<<E[2] << std::endl;
      // std::cout <<"H "<<H[0]<<", "<< H[1] <<", "<<H[2] << std::endl;
      // std::cout <<"vert "<<vert[0]<<", "<< vert[1] <<", "<<vert[2] << std::endl<<std::endl;
      //integral = integral + per_face_area_[i]*P;

      // Test Poynting vector integration
      std::vector<double> unit = { vert[0]/r, vert[1]/r, vert[2]/r};
      std::vector<double> P = (1/(2.0))
        *real(cross(E,vconj(H)));
      integral[0] = integral[0] + per_face_area_[i]*dot(P,unit);

    }
    return integral;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  
  std::vector<double> ShellGenerator::Integrate() {
    std::vector<double> integral = {0.0, 0.0, 0.0};
    //simple 
    for (unsigned int i=0; i<E_.size(); ++i) {
      auto E = E_[i];
      //auto H = 377.0*H_[i];
      auto H = H_[i];
      auto vert = vertices_[i];
      // Vector to unit product
      double r = norm(vert);
      std::vector<std::complex<double> > unit = { vert[0]/r, vert[1]/r, vert[2]/r};
      // std::cout << norm(unit) << std::endl;
      //const double pi = 3.1415926535897932384626433832795;
      std::vector<double> P = (1/(2.0))
        *real(
              dot(unit,E)*vconj(E) +
              dot(unit,H)*vconj(H) +
              (-1.0/2.0)*(dot(E,vconj(E))
                          +dot(H,vconj(H))
                          )*unit
              );
      // auto    std::cout <<"E "<<E[0]<<", "<< E[1] <<", "<<E[2] << std::endl;
      // std::cout <<"H "<<H[0]<<", "<< H[1] <<", "<<H[2] << std::endl;
      // std::cout <<"vert "<<vert[0]<<", "<< vert[1] <<", "<<vert[2] << std::endl<<std::endl;
      integral = integral + per_vertice_area_*P;
    }
    return integral;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector< std::vector<double> > ShellGenerator::GetVerticesT() {
    std::vector< std::vector<double> > vertices_t;
    vertices_t.resize(3); 
    for(const auto vert : vertices_){
      vertices_t[0].push_back(vert[0]);
      vertices_t[1].push_back(vert[1]);
      vertices_t[2].push_back(vert[2]);
    }
    return vertices_t;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector< std::vector<double> > ShellGenerator::GetFaceCentersT() {
    EvalFaces();
    std::vector< std::vector<double> > vertices_t;
    vertices_t.resize(3); 
    for(const auto vert : face_centers_){
      vertices_t[0].push_back(vert[0]);
      vertices_t[1].push_back(vert[1]);
      vertices_t[2].push_back(vert[2]);
    }
    return vertices_t;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double ShellGenerator::norm(std::vector<double> a){
    double norm_value = 0;
    for (auto coord:a)
      norm_value += pow2(coord);
    return std::sqrt(norm_value);
   }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
   void ShellGenerator::Rescale(double scale) {
     for(auto& vert : vertices_){
       double factor = norm(vert);
       //std::cout<< factor <<std::endl;
       for (auto &coord:vert) {
         coord = coord*scale/factor;
       }
       //std::cout << " " << norm(vert) << " ";
     }
     const double pi = 3.1415926535897932384626433832795;
     double area = 4.0*pi*pow2(scale); 
     //face_area_ = area/faces_.size();
     per_vertice_area_ = area/vertices_.size();
     //std::cout << "Per verice area: " << per_vertice_area_ << std::endl;
   }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
   void ShellGenerator::PrintVerts() {
   std::cout << "Verts coords:" << std::endl;
    for(auto vert : vertices_){
      std::cout <<"(";
      for (auto coord:vert) std::cout<<coord<<",";
      std::cout <<"),";      
    }
    std::cout << std::endl;
   }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::EvalFaces() {
    face_centers_.clear();
    per_face_area_.clear();
    for (auto face : faces_) {
      std::set<long unsigned int> edge_points;
      for (auto edge : face) {
        edge_points.insert(edges_[edge][0]);
        edge_points.insert(edges_[edge][1]);
      }
      std::vector<double> mid_point({0.0, 0.0, 0.0});
      for (auto point : edge_points) {
        mid_point = mid_point + vertices_[point];
      }
      mid_point = mid_point/3.0;
      face_centers_.push_back(mid_point);
      std::vector<long unsigned int> v_edge_points( edge_points.begin(),
                                                    edge_points.end() );
      auto vec_a = vertices_[v_edge_points[1]] - vertices_[v_edge_points[0]];
      auto vec_b = vertices_[v_edge_points[2]] - vertices_[v_edge_points[0]];
      auto area = norm(cross(vec_a, vec_b))/2.0;
      per_face_area_.push_back(area);
      //std::cout << "Area " <<area<<std::endl;
    }  // end for face in faces_
    double total_flat_area = 0.0;
    for (auto face:per_face_area_)
      total_flat_area += face;
    auto scale = norm(vertices_[0]);
    const double pi = 3.1415926535897932384626433832795;
    double area = 4.0*pi*pow2(scale); 
    face_area_ = area/faces_.size();
    double area_scale = area/total_flat_area;
    for (auto& face:per_face_area_)
      face *= area_scale;
    
    
    for(auto& vert : face_centers_){
       double factor = norm(vert);
       //std::cout<< factor <<std::endl;
       for (auto &coord:vert) {
         coord*=scale/factor;
       }
       //std::cout << " " << norm(vert) << " ";
     }

    // std::cout << "total face centers: " << face_centers_.size()
    //           << " scale " << scale
    //           << " face_norm " << norm(face_centers_[0])
    //           << " area-int " << face_area_
    //           << std::endl;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::Refine() {
    for (auto &edge : edges_) {
      auto p0 = vertices_[edge[0]];
      auto p1 = vertices_[edge[1]];
      std::vector<double> new_point = {
        (p0[0]+p1[0])/2.0,
        (p0[1]+p1[1])/2.0,
        (p0[2]+p1[2])/2.0};
      vertices_.push_back(new_point);
      edge.push_back(vertices_.size()-1);  // the last index is for the new mid-point
    }
    std::cout << "new verts: " << vertices_.size() <<std::endl;
    // std::cout << "extended edges:" <<std::endl;
    // for (auto edge : edges_) {
    //   std::cout<< "\t"<< edge[0]<< "\t"<< edge[1]<< "\t"<< edge[2]<<std::endl;
    // }
    refined_edges_.clear();
    for (auto edge : edges_) {
      // Important! New (refined) point goes the last!
      std::vector<long unsigned int> edge_a = {edge[0],edge[2]};
      std::vector<long unsigned int> edge_b = {edge[1],edge[2]};
      refined_edges_.push_back(edge_a);
      refined_edges_.push_back(edge_b);
    }
    // Now we need to count edges inside old faces.
    refined_faces_.clear();
    for (auto face : faces_) {
      auto point_a = edges_[face[0]][2];
      auto point_b = edges_[face[1]][2];
      auto point_c = edges_[face[2]][2];
      // std::cout << "\tedges_old: " <<face[0]<<" "<<face[1]<<" "<<face[2]<<" "<<std::endl;
      // std::cout << "\tpoints_old_edge0: " <<edges_[face[0]][0]<<" "<<edges_[face[0]][1]<<" "<<edges_[face[0]][2]<<" "<<std::endl;
      // std::cout << "\tpoints_old_edge0: " <<edges_[face[1]][0]<<" "<<edges_[face[1]][1]<<" "<<edges_[face[1]][2]<<" "<<std::endl;
      // std::cout << "\tpoints_old_edge0: " <<edges_[face[2]][0]<<" "<<edges_[face[2]][1]<<" "<<edges_[face[2]][2]<<" "<<std::endl;
      // std::cout<<"\trefined points: "<<point_a<<" "<<point_b<<" "<<point_c<<std::endl;
      std::vector<long unsigned int> edge_c = {point_a, point_b};
      std::vector<long unsigned int> edge_a = {point_b, point_c};
      std::vector<long unsigned int> edge_b = {point_c, point_a};
      refined_edges_.push_back(edge_a);
      auto edge_a_index = refined_edges_.size()-1;
      refined_edges_.push_back(edge_b);
      auto edge_b_index = refined_edges_.size()-1;
      refined_edges_.push_back(edge_c);
      auto edge_c_index = refined_edges_.size()-1;

      /*
      //                    /\      contrcloсkwise                         
      //                   c  1                                    
      //                  0    b                                   
      //                 /__a___\   edge_a                         
      //                /\      /\                                 
      //               c  \    /  0                                
      //              1    \  /    b                               
      //    edge_0a  /__0a__\/__1a__\ edge_1a                      
      //
      // remember! In edge_0a the refined point is [1], etc.
      */
      auto edge_0a = refined_edges_[2*face[0]];
      auto edge_1a = refined_edges_[2*face[0]+1];
      auto edge_0b = refined_edges_[2*face[1]];
      auto edge_1b = refined_edges_[2*face[1]+1];
      auto edge_0c = refined_edges_[2*face[2]];
      auto edge_1c = refined_edges_[2*face[2]+1];
      auto edge_0a_index = 2*face[0];
      auto edge_1a_index = 2*face[0]+1;
      auto edge_0b_index = 2*face[1];
      auto edge_1b_index = 2*face[1]+1;
      auto edge_0c_index = 2*face[2];
      auto edge_1c_index = 2*face[2]+1;
      // Orient:
      // Try contrcloсkwise:
      bool isClockwise = false, is_b_swapped = false, is_c_swapped=false;
      
      if (edge_0a[0]!=edge_1c[0]) {
        edge_1c.swap(edge_0c);
        is_c_swapped = !is_c_swapped;
      }
      if (edge_1a[0]!=edge_0b[0]) {
        edge_0b.swap(edge_1b);
        is_b_swapped = !is_b_swapped;
      }
      if (edge_1b[0]!=edge_0c[0]) {
        isClockwise = true;
        //Try clockwise:
        if (edge_0a[0]!=edge_1b[0]) {
          edge_1b.swap(edge_0b);
          is_b_swapped = !is_b_swapped;
        }
        if (edge_1a[0]!=edge_0c[0]) {
          edge_0c.swap(edge_1c);
          is_c_swapped = !is_c_swapped;
        }
        if (edge_1c[0]!=edge_0b[0])
          throw std::invalid_argument("Error! Unable to orient edges of refined face!\n");
      }
      if (is_b_swapped) {
       edge_1b_index = 2*face[1];
       edge_0b_index = 2*face[1]+1;
      }
      if (is_c_swapped) {
       edge_1c_index = 2*face[2];
       edge_0c_index = 2*face[2]+1;
      }

      /*
      //                    /\      clockwise                               
      //                   b  1                                    
      //                  0    c                                   
      //                 /__a___\   edge_a                         
      //                /\      /\                                 
      //               b  \    /  0                                
      //              1    \  /    c                               
      //    edge_0a  /__0a__\/__1a__\ edge_1a                      
      //
      */
      //Build new facets:
      // if isClockwise 
      std::vector<long unsigned int> face1({edge_0a_index, edge_1b_index, edge_c_index});
      std::vector<long unsigned int> face2({edge_1a_index, edge_0c_index, edge_b_index});
      std::vector<long unsigned int> face3({edge_0b_index, edge_1c_index, edge_a_index});
      std::vector<long unsigned int> face4({edge_a_index, edge_b_index, edge_c_index});
      if (!isClockwise) {
        face1 = std::vector<long unsigned int>({edge_0a_index, edge_1c_index, edge_b_index});
        face2 = std::vector<long unsigned int>({edge_1a_index, edge_0b_index, edge_c_index});
        face3 = std::vector<long unsigned int>({edge_1b_index, edge_0c_index, edge_a_index});
        face4 = std::vector<long unsigned int>({edge_a_index, edge_b_index, edge_c_index});
      }
      // std::cout<< "\tface1\t"<< face1[0]<< "\t"<< face1[1]<< "\t"<< face1[2]<<std::endl;
      // std::cout<< "\tface2\t"<< face2[0]<< "\t"<< face2[1]<< "\t"<< face2[2]<<std::endl;
      // std::cout<< "\tface3\t"<< face3[0]<< "\t"<< face3[1]<< "\t"<< face3[2]<<std::endl;
      // std::cout<< "\tface4\t"<< face4[0]<< "\t"<< face4[1]<< "\t"<< face4[2]<<std::endl;
      refined_faces_.push_back(face1);
      refined_faces_.push_back(face2);
      refined_faces_.push_back(face3);
      refined_faces_.push_back(face4);
      // std::cout<<"ref edges size: " << refined_edges_.size()<< std::endl;
      // std::cout << "Face1 points: "<< refined_edges_[face1[0]][0]
      //           << " " << refined_edges_[face1[0]][1] << "  --  "
      //           << refined_edges_[face1[1]][0]
      //           << " " << refined_edges_[face1[1]][1] << "  --  "
      // << refined_edges_[face1[2]][0]
      //                << " " << refined_edges_[face1[2]][1] << std::endl;

    } // end for faces_
    
    std::cout << "new edges: " << refined_edges_.size() <<std::endl;
    std::cout << "new faces: " << refined_faces_.size() <<std::endl;
    edges_.clear();
    edges_ = refined_edges_;
    // std::cout << "edges:" <<std::endl;
    // for (auto edge : edges_) {
    //   std::cout<< " "<< edge[0]<< "\t"<< edge[1]<<std::endl;
    // }
    faces_.clear();
    faces_ = refined_faces_;
    
    //Rescale(1.0);
    //GenerateEdges();
    // GenerateFaces();
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::GenerateFacesInit() {
    faces_.clear();
    for (unsigned int i = 0; i < edges_.size(); ++i) {
      const auto ie = edges_[i];
      for(unsigned int j = i + 1; j < edges_.size(); ++j) {
        const auto je = edges_[j];
        for(unsigned int k = j + 1; k < edges_.size(); ++k) {
          const auto ke = edges_[k];
          std::set<long unsigned int> edges = {ie[0],ie[1],
                                 je[0],je[1],
                                 ke[0],ke[1]};
          if (edges.size() != 3) continue;
          std::vector<long unsigned int> face({i,j,k});
          // std::cout << ie[0]<<"-"<<ie[1] << ":"
          //             << je[0]<<"-"<<je[1] << ":"
          //             << ke[0]<<"-"<<ke[1] 
          //             << std::endl;
          // std::cout << face[0]<<"-"<<face[1] << "-"<<face[2]<<std::endl;
          faces_.push_back(face);
        }
      }
    }
    std::cout << "Faces: "<< faces_.size() <<std::endl;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::GenerateEdgesInit() {
    std::cout << "Vertices: "<< vertices_.size() <<std::endl;
    edges_.clear();
    EvalEdgeLength();
    for (unsigned int i = 0; i < vertices_.size(); ++i)
      for(unsigned int j = i + 1; j < vertices_.size(); ++j) {
        //std::cout << i<< " " << j<<  " == "<< dist(vertices_[i],vertices_[j]) <<std::endl;
        if (dist(vertices_[i],vertices_[j]) > 1.000001*edge_length_) continue;
        std::vector<long unsigned int> edge = {i,j};
        // std::cout << i<< " " << j << " : i=(";
        // for (auto v : vertices_[i]) std::cout << v <<",";
        // std::cout<<")  j=(";
        // for (auto v : vertices_[j]) std::cout << v <<",";
        // std::cout<<")"<<std::endl;
        edges_.push_back(edge);
      }
    std::cout << "Edges: "<< edges_.size() <<std::endl;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::EvalEdgeLength() {
    auto p0 = vertices_[0];
    std::vector<double> zero(p0.size(), 0.0);
    double min_dist = 42.0*dist(zero, p0);
    for (auto point : vertices_) {
      if (point == p0) continue;
      double new_dist = dist(p0, point);
      if (new_dist < min_dist) min_dist = new_dist;
    }
    //std::cout << "Edge length = " << min_dist << std::endl;
    edge_length_ = min_dist;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double ShellGenerator::dist(std::vector<double> a, std::vector<double> b) {
    unsigned int len = a.size();
    if (b.size() != len)
      throw std::invalid_argument("Error! Vector need to be the same size!\n");
    double distance = 0.0;
    for (unsigned int i = 0; i<len; ++i) {
      distance += pow2(a[i]-b[i]);
    }
    return std::sqrt(distance);
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // @brief set up regular icosahedron
  void ShellGenerator::SetInitialVerticesIcosahedron() {
    double a = 0.0;
    double b = 1.0;
    double c = (1+std::sqrt(5.0))/2.0;
    std::vector< std::vector<double> > points = {
      {a, b, c},
      {a, b,-c},
      {a,-b, c},
      {a,-b,-c},
      
      { b, c,a},
      { b,-c,a},
      {-b, c,a},
      {-b,-c,a},

      { c,a, b},
      {-c,a, b},
      { c,a,-b},
      {-c,a,-b}
    };
    vertices_ = std::move(points);
    //Rescale(1.0);

    //std::vector< std::vector<double> > points_debug = {{1,0,0},{-1,0,0}};
    //std::vector< std::vector<double> > points_debug = {{0,1,0},{0,-1,0}};
    //std::vector< std::vector<double> > points_debug = {{1,1,0},{1,-1,0},{-1,-1,0},{-1,1,0}};
    //std::vector< std::vector<double> > points_debug = {{0,1,1},{0,1,-1},{0,-1,-1},{0,-1,1}};
    //std::vector< std::vector<double> > points_debug = {};
    // std::vector< std::vector<double> > points_debug = {{0,0,1},{0,0,-1}};

    //vertices_ = std::move(points_debug);

    // for (auto v : vertices_) {
    //   for (auto p : v)
    //     std::cout<< p<<std::endl;
    // }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // @brief set up regular icosahedron
  void ShellGenerator::SetInitialVerticesTetrahedron() {
    double a = 1.0;
    double b = 1.0;
    double c = 1.0;
    std::vector< std::vector<double> > points = {
      {a, b, c},
      {a, -b,-c},
      {-a,-b, c},
      {-a, b,-c}      
    };
    vertices_ = std::move(points);
    //Rescale(1.0);


    //vertices_ = std::move(points_debug);

    // for (auto v : vertices_) {
    //   for (auto p : v)
    //     std::cout<< p<<std::endl;
    // }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::Init() {
    SetInitialVerticesIcosahedron();
    //SetInitialVerticesTetrahedron();
    GenerateEdgesInit();
    GenerateFacesInit();
  }

}  // end of namespace read_spectra

