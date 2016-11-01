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
  // template<typename T, typename A> inline
  // std::vector<T> cross(std::vector<T,A> const& a,
  //       std::vector<T,A> const& b) {
  //   std::vector<T> c = {
  //     a[1]*b[2]-a[2]*b[1],
  //     a[2]*b[0]-a[0]*b[2],
  //     a[0]*b[1]-a[1]*b[0]
  //   };
  //   return c;
  // }  
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
  double ShellGenerator::IntegrateGaussSimple(double charge, double shift) {
    double integral = 0.0;
    // return field at point p from the charge, located at (shift, 0,0)
    auto field = [](double charge, double shift, std::vector<double> p){
      double r = std::sqrt(pow2(p[0]-shift) + pow2(p[1]) + pow2(p[2]) );
      const double pi = 3.1415926535897932384626433832795;
      double ampl = charge/(4.0*pi*pow2(r));      
      std::vector<double> field = {ampl*(p[0]-shift)/r, ampl*(p[1])/r, ampl*(p[2])/r};
      return field;
    };
    //simple 
    for (auto vert :vertices_) {
      auto E0 = field(charge, shift, vert);
      // Vector to unit product
      double r = norm(vert);
      std::vector<double> unit = { vert[0]/r, vert[1]/r, vert[2]/r};
      // std::cout << norm(unit) << std::endl;
      for (int i =0; i < 3; ++i) 
        integral += per_vertice_area_*unit[i]*E0[i];
    }
    return integral;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double ShellGenerator::IntegrateGauss(double charge, double shift) {
    if (faces_.size() == 0)
      throw std::invalid_argument("Error! Faces were not initialized!\n");
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
  
  std::vector<double> ShellGenerator::Integrate() {
    std::vector<double> integral = {0.0, 0.0, 0.0};
    //simple 
    for (unsigned int i=0; i<E_.size(); ++i) {
      auto E = E_[i];
      //auto H = 377.0*H_[i];
      auto H = H_[i];
      auto vert = vertices_[i];
      std::cout <<"E "<<E[0]<<", "<< E[1] <<", "<<E[2] << std::endl;
      std::cout <<"H "<<H[0]<<", "<< H[1] <<", "<<H[2] << std::endl;
      std::cout <<"vert "<<vert[0]<<", "<< vert[1] <<", "<<vert[2] << std::endl<<std::endl;
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
       double factor = scale/norm(vert);
       //std::cout<< factor <<std::endl;
       for (auto &coord:vert) {
         coord*=factor;
       }
     }
     const double pi = 3.1415926535897932384626433832795;
     double area = 4.0*pi*pow2(scale); 
     face_area_ = area/faces_.size();
     per_vertice_area_ = area/vertices_.size(); 
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
  void ShellGenerator::Refine() {
    for (auto edge : edges_) {
      auto p0 = vertices_[edge[0]];
      auto p1 = vertices_[edge[1]];
      std::vector<double> new_point = {
        (p0[0]+p1[0])/2.0,
        (p0[1]+p1[1])/2.0,
        (p0[2]+p1[2])/2.0};
      vertices_.push_back(new_point);      
    }
    GenerateEdges();
    // GenerateFaces();
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::GenerateFaces() {
    faces_.clear();
    for (unsigned int i = 0; i < edges_.size(); ++i) {
      const auto ie = edges_[i];
      for(unsigned int j = i + 1; j < edges_.size(); ++j) {
        const auto je = edges_[j];
        for(unsigned int k = j + 1; k < edges_.size(); ++k) {
          const auto ke = edges_[k];
          std::set<long int> edges = {ie[0],ie[1],
                                 je[0],je[1],
                                 ke[0],ke[1]};
          if (edges.size() != 3) continue;
          std::vector<long int> face(edges.begin(), edges.end());
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
  void ShellGenerator::GenerateEdges() {
    std::cout << "Vertices: "<< vertices_.size() <<std::endl;
    edges_.clear();
    EvalEdgeLength();
    for (unsigned int i = 0; i < vertices_.size(); ++i)
      for(unsigned int j = i + 1; j < vertices_.size(); ++j) {
        if (dist(vertices_[i],vertices_[j]) > 1.001*edge_length_) continue;
        std::vector<long int> edge = {i,j};
        // std::cout << i<< " " << j<<std::endl;
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
    std::cout << "Edge length = " << min_dist << std::endl;
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
  void ShellGenerator::SetInitialVertices() {
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

    //std::vector< std::vector<double> > points_debug = {{1,0,0},{-1,0,0}};
    //std::vector< std::vector<double> > points_debug = {{0,1,0},{0,-1,0}};
    std::vector< std::vector<double> > points_debug = {{1,1,0},{-1,-1,0},{-1,1,0},{1,-1,0}};
    //std::vector< std::vector<double> > points_debug = {};
    // std::vector< std::vector<double> > points_debug = {{0,0,1},{0,0,-1}};
    vertices_ = std::move(points_debug);

    // for (auto v : vertices_) {
    //   for (auto p : v)
    //     std::cout<< p<<std::endl;
    // }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::Init() {
    SetInitialVertices();
    GenerateEdges();
    //GenerateFaces();
  }

}  // end of namespace read_spectra

