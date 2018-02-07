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
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "shell-generator.hpp"
namespace shell_generator {
  struct KahanAccumulation
  {
    double sum;
    double correction;
  };
 
  KahanAccumulation KahanSum(KahanAccumulation accumulation, double value)
  {
    KahanAccumulation result;
    double y = value - accumulation.correction;
    double t = accumulation.sum + y;
    result.correction = (t - accumulation.sum) - y;
    result.sum = t;
    return result;
  }
  
  double ShellGenerator::norm(std::vector<double> a){
    double norm_value = 0;
    for (auto coord:a)
      norm_value += pow2(coord);
    return std::sqrt(norm_value);
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
  std::vector<double> ShellGenerator::GetDipoleField(std::vector<double> point) {
    double charge = 3.14;
    double dist = 0;
    std::vector< std::complex<double> > zero (3,std::complex<double>(0.0,0.0));
    auto dim = vertices_.size();
    E_.clear();
    H_.clear();
    E_ = std::vector< std::vector< std::complex<double> > > (dim, zero);
    H_ = std::vector< std::vector< std::complex<double> > > (dim, zero);

    std::vector<double> shift = {dist, 0.0, 0.0};
    for (long unsigned int i=0; i< dim; ++i) {
      auto vert = vertices_[i];
      double r = norm(vert);
      double r2 = norm(vert-shift);
      std::vector<std::complex<double> > unit = { vert[0]/r, vert[1]/r, vert[2]/r};
      std::vector<std::complex<double> > unit2 = { (vert[0]-dist)/r2, vert[1]/r2, vert[2]/r2};
    const double pi = 3.1415926535897932384626433832795;
  //const double pi = 3.1415926535897932384626433832795;
      double ampl = charge/(4.0*pi*pow2(r));      
      double ampl2 = charge/(4.0*pi*pow2(r2));      
      E_[i] = ampl*unit + ampl2*unit2;
      // for (auto E:E_[i]) std::cout << E << " ";
      // std::cout<< "  <-  "<<norm(real(E_[i]))<<std::endl;
      
    }
      return point;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  
  std::vector<double> ShellGenerator::IntegrateByFacesQuadrature2() {
    std::vector<double> integral = {0.0, 0.0, 0.0};
    std::vector<double> int1, int2, int3;
    for (auto face : faces_) {
      std::vector< std::vector<double> > quadrature_points;
      for (auto edge : face) {
        auto point_a = vertices_[edges_[edge][0]];
        auto point_b = vertices_[edges_[edge][1]];
        auto mid_point = 1/2.0*(point_a+point_b);
        // Important! Rescale quadrature points  to the exact position
        // on integration sphere
        const double factor = norm(mid_point);
        for (auto &coord:mid_point) coord = coord*scale_/factor;
        quadrature_points.push_back(mid_point);
      }
      for (auto point:quadrature_points){
        auto value = (*ValueAtPoint)(point);
        //integral = integral + (1/3.0)*face_area_*(*ValueAtPoint)(point);
        auto pre_sum = (1/3.0)*face_area_*value;
        //std::cout << face_area_ << " : "<<value[0]<<std::endl;    

        int1.push_back(pre_sum[0]);
        int2.push_back(pre_sum[1]);
        int3.push_back(pre_sum[2]);
      }
    }
    { KahanAccumulation init = {0};
    KahanAccumulation result =
        std::accumulate(int1.begin(), int1.end(), init, KahanSum);
    integral[0] = result.sum; }
    { KahanAccumulation init = {0};
    KahanAccumulation result =
        std::accumulate(int2.begin(), int2.end(), init, KahanSum);
    integral[1] = result.sum; }
    { KahanAccumulation init = {0};
    KahanAccumulation result =
        std::accumulate(int3.begin(), int3.end(), init, KahanSum);
    integral[2] = result.sum; }
    
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
   void ShellGenerator::Rescale(double scale) {
     scale_ = scale;
     for(auto& vert : vertices_){
       double factor = norm(vert);
       //std::cout<< factor <<std::endl;
       for (auto &coord:vert) {
         coord = coord*scale/factor;
       }
       //std::cout << " " << norm(vert) << " ";
     }
   const double pi = 3.1415926535897932384626433832795;
 //const double pi = 3.1415926535897932384626433832795;
     double area = 4.0*pi*pow2(scale); 
     face_area_ = area/faces_.size();
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
  //const double pi = 3.1415926535897932384626433832795;
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
      // the last index is for the new mid-point
      vertices_.push_back(new_point);
      // now it will be a new point on the edge, the last (numbered [2])
      // entry in the egde points list.
      edge.push_back(vertices_.size()-1);  
    }
    //std::cout << "new verts: " << vertices_.size() <<std::endl;
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
    
    // std::cout << "new edges: " << refined_edges_.size() <<std::endl;
    // std::cout << "new faces: " << refined_faces_.size() <<std::endl;
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
    //std::cout << "Faces: "<< faces_.size() <<std::endl;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ShellGenerator::GenerateEdgesInit() {
    //std::cout << "Vertices: "<< vertices_.size() <<std::endl;
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
    //std::cout << "Edges: "<< edges_.size() <<std::endl;
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
    scale_ = std::sqrt( pow2(a) + pow2(b) + pow2(c));
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

