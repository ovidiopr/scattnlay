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
  double ShellGenerator::Integrate() {
    double integral = 0.0;
    auto E = E_[0];
    auto H = H_[0];
    for(const auto face : faces_){
      for (int i = 0; i<3; ++i) {
        auto mean = (E_[face[0]][i]+E_[face[1]][i]+E_[face[2]][i])/3.0;
        E[i] = mean;
        mean = (H_[face[0]][i]+H_[face[1]][i]+H_[face[2]][i])/3.0;
        H[i] = mean;
      }
      integral += // std::abs(E[0]*H[0])*
        face_area_;
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
    for (int i = 0; i < edges_.size(); ++i) {
      const auto ie = edges_[i];
      for(int j = i + 1; j < edges_.size(); ++j) {
        const auto je = edges_[j];
        for(int k = j + 1; k < edges_.size(); ++k) {
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
    for (int i = 0; i < vertices_.size(); ++i)
      for(int j = i + 1; j < vertices_.size(); ++j) {
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
    for (int i = 0; i<len; ++i) {
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

