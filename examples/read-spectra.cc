/**
 * @file   read-spectra.cc
 * @author Konstantin Ladutenko <kostyfisik at gmail (.) com>
 * @date   Wed Mar 11 11:51:26 2015
 * 
 * @copyright 2015 Konstantin Ladutenko
 *
 * @brief  Read complex spectra from file in format 'WL real imag'
 * 
 * read-spectra is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * read-spectra is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with read-spectra.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include <algorithm>
#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <iostream>
#include <vector>
#include "read-spectra.h"
namespace read_spectra {
  template<class T> inline T pow2(const T value) {return value*value;}
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  ReadSpectra& ReadSpectra::ReadFromFile(std::string fname) {
    //std::cout<<"Reading file: "<< fname << std::endl;
    std::ifstream infile(fname.c_str());
    data_.clear();
    std::string line;
    while (std::getline(infile, line))
      {
	if (line.front() == '#') continue; //do not read comments
	if (line.find('#') != std::string::npos) 
	  throw std::invalid_argument("Error! Comments should be marked with # in the begining of the line!\n");
	std::istringstream iss(line);	
	double wl, re, im;
	if (!(iss >> wl >> re >> im)) throw std::invalid_argument("Error! Unexpected format of the line!\n");
	data_.push_back(std::vector<double>({wl,re,im}));
	//std::cout<<wl<<' '<<re<<' '<<im<<std::endl;
      }  // end of wile reading file 
    std::sort(data_.begin(), data_.end(),
	      [](const std::vector<double>& a, const std::vector<double>& b) {
		return a.front() < b.front();
	      });
    return *this;
  }  // end of void ReadSpectra::ReadFromFile(std::string fname)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// Cut the spectra to the range and convert it to std::complex<double>
  ReadSpectra& ReadSpectra::ResizeToComplex(double from_wl, double to_wl, int samples) {
    if (data_.size() < 2) throw std::invalid_argument("Nothing to resize!/n");
    if (data_.front()[0] > from_wl || data_.front()[0] > to_wl ||
	data_.back()[0] < from_wl || data_.back()[0] < to_wl ||
	from_wl > to_wl)
      throw std::invalid_argument("Invalid range to resize spectra!/n");
    if (samples < 1) throw std::invalid_argument("Not enough samples!/n");
    std::vector<double> wl_sampled(samples, 0.0);
    if (samples == 1) {
      wl_sampled[0] = (from_wl + to_wl)/2.0;
    } else {
      for (int i =0; i<samples; ++i)
	wl_sampled[i] = from_wl
	  + (to_wl-from_wl)*static_cast<double>(i)/static_cast<double>(samples-1);
    }  // end of setting wl_sampled
    data_complex_.clear();
    int j = 0;
    for (int i = 0; i < data_.size(); ++i) {
      const double& wl_i = data_[i][0];
      const double& wl_s = wl_sampled[j];
      if (wl_i < wl_s) continue;
      else {
	const double& wl_prev = data_[i-1][0];	
	const double& re_prev = data_[i-1][1];
	const double& im_prev = data_[i-1][2];
	const double& re_i = data_[i][1];
	const double& im_i = data_[i][2];
	// Linear approximation
	double re_s = re_i - (re_i-re_prev)*(wl_i-wl_s)/(wl_i-wl_prev);
	double im_s = im_i - (im_i-im_prev)*(wl_i-wl_s)/(wl_i-wl_prev);

	auto tmp = std::make_pair(wl_s, std::complex<double>(re_s,im_s));
	data_complex_.push_back(tmp);	
	       
	++j;
	--i; // Next sampled point(j) can be in the same i .. i-1 region
	// All sampled wavelengths has a value
	if (j >= wl_sampled.size()) break;  
      }
    }
    if (data_complex_.size() == 0)
      throw std::invalid_argument("No points in spectra for the defined range!/n");
    if (data_complex_.size() != samples)
      throw std::invalid_argument("Was not able to get all samples!/n");
    return *this;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// from relative permittivity to refractive index
  ReadSpectra& ReadSpectra::ToIndex() {
    data_complex_index_.clear();
    for (auto row : data_complex_) {
      const double wl = row.first;
      const double e1 = row.second.real();
      const double e2 = row.second.imag();
      const double n = std::sqrt( (std::sqrt(pow2(e1)+pow2(e2)) + e1) /2.0 );
      const double k = std::sqrt( (std::sqrt(pow2(e1)+pow2(e2)) - e1) /2.0 );
      auto tmp = std::make_pair(wl, std::complex<double>(n,k));
      data_complex_index_.push_back(tmp);	
    }
    return *this;
  }

  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void ReadSpectra::PrintData() {
    if (data_complex_.size() == 0) 
      throw std::invalid_argument("Nothing to print!");
    for (auto row : data_complex_) {
      printf("wl:%g\tre:%g\tim:%g\n", row.first, row.second.real(),
	     row.second.imag());
    }  // end of for each row
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //

}  // end of namespace read_spectra

