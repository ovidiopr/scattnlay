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
  double ShellGenerator::IntegrateGauss(double charge, double shift) {
    double integral = 0.0;
    auto E = E_[0];
    auto H = H_[0];

    // return field at point p from the charge, located at (shift, 0,0)
    auto field = [](double charge, double shift, std::vector<double> p){
      double r = std::sqrt(pow2(p[0]-shift) + pow2(p[1]) + pow2(p[2]) );
      double ampl = charge/pow2(r);
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
      std::cout << norm(unit) << std::endl;
      for (int i =0; i < 3; ++i) 
        integral += face_area_*unit[i]*mean_vector[i];
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
     face_area_ = 4.0*pi*pow2(scale)/faces_.size();
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
    GenerateFaces();
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



  // ShellGenerator& ShellGenerator::ReadFromFile(std::string fname) {
  //   //std::cout<<"Reading file: "<< fname << std::endl;
  //   std::ifstream infile(fname.c_str());
  //   data_.clear();
  //   std::string line;
  //   while (std::getline(infile, line))
  //     {
  //       if (line.front() == '#') continue; //do not read comments
  //       if (line.find('#') != std::string::npos) 
  //         throw std::invalid_argument("Error! Comments should be marked with # in the begining of the line!\n");
  //       std::istringstream iss(line);	
  //       double wl, re, im;
  //       if (!(iss >> wl >> re >> im)) throw std::invalid_argument("Error! Unexpected format of the line!\n");
  //       data_.push_back(std::vector<double>({wl,re,im}));
  //       //std::cout<<wl<<' '<<re<<' '<<im<<std::endl;
  //     }  // end of wile reading file 
  //   std::sort(data_.begin(), data_.end(),
  //             [](const std::vector<double>& a, const std::vector<double>& b) {
  //       	return a.front() < b.front();
  //             });
  //   return *this;
  // }  // end of void ShellGenerator::ReadFromFile(std::string fname)
  // // ********************************************************************** //
  // // ********************************************************************** //
  // // ********************************************************************** //
  // /// Cut the spectra to the range and convert it to std::complex<double>
  // ShellGenerator& ShellGenerator::ResizeToComplex(double from_wl, double to_wl, int samples) {
  //   if (data_.size() < 2) throw std::invalid_argument("Nothing to resize!/n");
  //   if (data_.front()[0] > from_wl || data_.front()[0] > to_wl ||
  //       data_.back()[0] < from_wl || data_.back()[0] < to_wl ||
  //       from_wl > to_wl)
  //     throw std::invalid_argument("Invalid range to resize spectra!/n");
  //   if (samples < 1) throw std::invalid_argument("Not enough samples!/n");
  //   std::vector<double> wl_sampled(samples, 0.0);
  //   if (samples == 1) {
  //     wl_sampled[0] = (from_wl + to_wl)/2.0;
  //   } else {
  //     for (int i =0; i<samples; ++i)
  //       wl_sampled[i] = from_wl
  //         + (to_wl-from_wl)*static_cast<double>(i)/static_cast<double>(samples-1);
  //   }  // end of setting wl_sampled
  //   data_complex_.clear();
  //   int j = 0;
  //   for (int i = 0; i < data_.size(); ++i) {
  //     const double& wl_i = data_[i][0];
  //     const double& wl_s = wl_sampled[j];
  //     if (wl_i < wl_s) continue;
  //     else {
  //       const double& wl_prev = data_[i-1][0];	
  //       const double& re_prev = data_[i-1][1];
  //       const double& im_prev = data_[i-1][2];
  //       const double& re_i = data_[i][1];
  //       const double& im_i = data_[i][2];
  //       // Linear approximation
  //       double re_s = re_i - (re_i-re_prev)*(wl_i-wl_s)/(wl_i-wl_prev);
  //       double im_s = im_i - (im_i-im_prev)*(wl_i-wl_s)/(wl_i-wl_prev);

  //       auto tmp = std::make_pair(wl_s, std::complex<double>(re_s,im_s));
  //       data_complex_.push_back(tmp);	
	       
  //       ++j;
  //       --i; // Next sampled point(j) can be in the same i .. i-1 region
  //       // All sampled wavelengths has a value
  //       if (j >= wl_sampled.size()) break;  
  //     }
  //   }
  //   if (data_complex_.size() == 0)
  //     throw std::invalid_argument("No points in spectra for the defined range!/n");
  //   if (data_complex_.size() != samples)
  //     throw std::invalid_argument("Was not able to get all samples!/n");
  //   return *this;
  // }
  // // ********************************************************************** //
  // // ********************************************************************** //
  // // ********************************************************************** //
  // /// from relative permittivity to refractive index
  // ShellGenerator& ShellGenerator::ToIndex() {
  //   data_complex_index_.clear();
  //   for (auto row : data_complex_) {
  //     const double wl = row.first;
  //     const double e1 = row.second.real();
  //     const double e2 = row.second.imag();
  //     const double n = std::sqrt( (std::sqrt(pow2(e1)+pow2(e2)) + e1) /2.0 );
  //     const double k = std::sqrt( (std::sqrt(pow2(e1)+pow2(e2)) - e1) /2.0 );
  //     auto tmp = std::make_pair(wl, std::complex<double>(n,k));
  //     data_complex_index_.push_back(tmp);	
  //   }
  //   return *this;
  // }

  // // ********************************************************************** //
  // // ********************************************************************** //
  // // ********************************************************************** //
  // void ShellGenerator::PrintData() {
  //   if (data_complex_.size() == 0) 
  //     throw std::invalid_argument("Nothing to print!");
  //   for (auto row : data_complex_) {
  //     printf("wl:%g\tre:%g\tim:%g\n", row.first, row.second.real(),
  //            row.second.imag());
  //   }  // end of for each row
  // }
  // // ********************************************************************** //
  // // ********************************************************************** //
  // // ********************************************************************** //

  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::Init() {
    SetInitialVertices();
    GenerateEdges();
    GenerateFaces();
  }

}  // end of namespace read_spectra

