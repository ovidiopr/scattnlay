///
/// @file   nmie.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Tue Sep  3 00:38:27 2013
/// @copyright 2013,2014,2015 Ladutenko Konstantin
///
/// nmie is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// nmie-wrapper is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with nmie-wrapper.  If not, see <http://www.gnu.org/licenses/>.
///
/// nmie uses nmie.c from scattnlay by Ovidio Pena
/// <ovidio@bytesfall.com> . He has an additional condition to 
/// his library:
//    The only additional condition is that we expect that all publications         //
//    describing  work using this software , or all commercial products             //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
///
/// @brief  Wrapper class around nMie function for ease of use
///
#include "nmie-wrapper.h"
#include <array>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace nmie {  
  //helpers
  template<class T> inline T pow2(const T value) {return value*value;}
  //#define round(x) ((x) >= 0 ? (int)((x) + 0.5):(int)((x) - 0.5))
  int round(double x) {
    return x >= 0 ? (int)(x + 0.5):(int)(x - 0.5);
  }  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  //emulate C call.
  int nMie_wrapper(int L, const std::vector<double>& x, const std::vector<std::complex<double> >& m, int nTheta, const std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
    
    if (x.size() != L || m.size() != L)
        throw std::invalid_argument("Declared number of layers do not fit x and m!");
    if (Theta.size() != nTheta)
        throw std::invalid_argument("Declared number of sample for Theta is not correct!");
    try {
      MultiLayerMie multi_layer_mie;  
      multi_layer_mie.SetWidthSP(x);
      multi_layer_mie.SetIndexSP(m);
      multi_layer_mie.SetAngles(Theta);
    
      multi_layer_mie.RunMieCalculations();
      
      *Qext = multi_layer_mie.GetQext();
      *Qsca = multi_layer_mie.GetQsca();
      *Qabs = multi_layer_mie.GetQabs();
      *Qbk = multi_layer_mie.GetQbk();
      *Qpr = multi_layer_mie.GetQpr();
      *g = multi_layer_mie.GetAsymmetryFactor();
      *Albedo = multi_layer_mie.GetAlbedo();
      S1 = multi_layer_mie.GetS1();
      S2 = multi_layer_mie.GetS2();
      //multi_layer_mie.GetFailed();
    } catch(const std::invalid_argument& ia) {
      // Will catch if  multi_layer_mie fails or other errors.
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      throw std::invalid_argument(ia);
      return -1;
    }  

    return 0;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int nField(const int L, const int pl, const std::vector<double>& x, const std::vector<std::complex<double> >& m, const int nmax, const int ncoord, const std::vector<double>& Xp_vec, const std::vector<double>& Yp_vec, const std::vector<double>& Zp_vec,  std::vector<std::vector<std::complex<double> > >& E, std::vector<std::vector<std::complex<double> > >& H) {
    if (x.size() != L || m.size() != L)
      throw std::invalid_argument("Declared number of layers do not fit x and m!");
    if (Xp_vec.size() != ncoord || Yp_vec.size() != ncoord || Zp_vec.size() != ncoord
        || E.size() != ncoord || H.size() != ncoord)
      throw std::invalid_argument("Declared number of coords do not fit Xp, Yp, Zp, E, or H!");
    for (auto f:E)
      if (f.size() != 3)
        throw std::invalid_argument("Field E is not 3D!");
    for (auto f:H)
      if (f.size() != 3)
        throw std::invalid_argument("Field H is not 3D!");
    try {
      MultiLayerMie multi_layer_mie;  
      //multi_layer_mie.SetPEC(pl);
      multi_layer_mie.SetWidthSP(x);
      multi_layer_mie.SetIndexSP(m);      
      multi_layer_mie.SetFieldPointsSP({Xp_vec, Yp_vec, Zp_vec});
      multi_layer_mie.RunFieldCalculations();
      E = multi_layer_mie.GetFieldE();
      H = multi_layer_mie.GetFieldH();
      //multi_layer_mie.GetFailed();
    } catch(const std::invalid_argument& ia) {
      // Will catch if  multi_layer_mie fails or other errors.
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      throw std::invalid_argument(ia);
      return - 1;
    }  

    return 0;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::GetFailed() {
    double faild_x = 9.42477796076938;
    //double faild_x = 9.42477796076937;
    std::complex<double> z(faild_x, 0.0);
    std::vector<int> nmax_local_array = {20, 100, 500, 2500};
    for (auto nmax_local : nmax_local_array) {
      std::vector<std::complex<double> > D1_failed(nmax_local + 1);
      // Downward recurrence for D1 - equations (16a) and (16b)
      D1_failed[nmax_local] = std::complex<double>(0.0, 0.0);
      const std::complex<double> zinv = std::complex<double>(1.0, 0.0)/z;
      for (int n = nmax_local; n > 0; n--) {
        D1_failed[n - 1] = double(n)*zinv - 1.0/(D1_failed[n] + double(n)*zinv);
      }
      printf("Faild D1[0] from reccurence (z = %16.14f, nmax = %d): %g\n",
             faild_x, nmax_local, D1_failed[0].real());
    }
    printf("Faild D1[0] from continued fraction (z = %16.14f): %g\n", faild_x,
           calcD1confra(0,z).real());
    //D1[nmax_] = calcD1confra(nmax_, z);
  
    
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMie::GetQext() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qext_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMie::GetQabs() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qabs_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> MultiLayerMie::GetQabs_channel() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qabs_ch_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> MultiLayerMie::GetQabs_channel_normalized() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    // std::vector<double> NACS(nmax_ - 1, 0.0);
    // double x2 = pow2(size_parameter_.back());
    // for (int i = 0; i < nmax_ - 1; ++i) {
    //   const int n = i + 1;
    //   NACS[i] = Qabs_ch_[i]*x2/(2.0*(2.0*static_cast<double>(n) + 1));
    //   // if (NACS[i] > 0.250000001)
    //   //         throw std::invalid_argument("Unexpected normalized absorption cross-section value!");
    // }
    //return NACS;    
    return Qabs_ch_norm_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMie::GetQsca() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qsca_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> MultiLayerMie::GetQsca_channel() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qsca_ch_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> MultiLayerMie::GetQsca_channel_normalized() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    // std::vector<double> NACS(nmax_ - 1, 0.0);
    // double x2 = pow2(size_parameter_.back());
    // for (int i = 0; i < nmax_ - 1; ++i) {
    //   const int n = i + 1;
    //   NACS[i] = Qsca_ch_[i]*x2/(2.0*(2.0*static_cast<double>(n) + 1.0));
    // }
    // return NACS;    
     return Qsca_ch_norm_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMie::GetQbk() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qbk_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMie::GetQpr() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qpr_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMie::GetAsymmetryFactor() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return asymmetry_factor_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMie::GetAlbedo() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return albedo_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<std::complex<double> > MultiLayerMie::GetS1() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return S1_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<std::complex<double> > MultiLayerMie::GetS2() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return S2_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::AddTargetLayer(double width, std::complex<double> layer_index) {
    isMieCalculated_ = false;
    if (width <= 0)
      throw std::invalid_argument("Layer width should be positive!");
    target_width_.push_back(width);
    target_index_.push_back(layer_index);
  }  // end of void  MultiLayerMie::AddTargetLayer(...)  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetTargetPEC(double radius) {
    isMieCalculated_ = false;
    if (target_width_.size() != 0 || target_index_.size() != 0)
      throw std::invalid_argument("Error! Define PEC target radius before any other layers!");
    // Add layer of any index...
    AddTargetLayer(radius, std::complex<double>(0.0, 0.0));
    // ... and mark it as PEC
    SetPEC(0.0);
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetCoatingIndex(std::vector<std::complex<double> > index) {
    isMieCalculated_ = false;
    coating_index_.clear();
    for (auto value : index) coating_index_.push_back(value);
  }  // end of void MultiLayerMie::SetCoatingIndex(std::vector<complex> index);  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetAngles(const std::vector<double>& angles) {
    isMieCalculated_ = false;
    theta_ = angles;
    // theta_.clear();
    // for (auto value : angles) theta_.push_back(value);
  }  // end of SetAngles()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetCoatingWidth(std::vector<double> width) {
    isMieCalculated_ = false;
    coating_width_.clear();
    for (auto w : width)
      if (w <= 0)
        throw std::invalid_argument("Coating width should be positive!");
      else coating_width_.push_back(w);
  }
  // end of void MultiLayerMie::SetCoatingWidth(...);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetWidthSP(const std::vector<double>& size_parameter) {
    isMieCalculated_ = false;
    size_parameter_.clear();
    double prev_size_parameter = 0.0;
    for (auto layer_size_parameter : size_parameter) {
      if (layer_size_parameter <= 0.0)
        throw std::invalid_argument("Size parameter should be positive!");
      if (prev_size_parameter > layer_size_parameter) 
        throw std::invalid_argument
          ("Size parameter for next layer should be larger than the previous one!");
      prev_size_parameter = layer_size_parameter;
      size_parameter_.push_back(layer_size_parameter);
    }
  }
  // end of void MultiLayerMie::SetWidthSP(...);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetIndexSP(const std::vector< std::complex<double> >& index) {
    isMieCalculated_ = false;
    //index_.clear();
    index_ = index;
    // for (auto value : index) index_.push_back(value);
  }  // end of void MultiLayerMie::SetIndexSP(...);  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetFieldPointsSP(const std::vector< std::vector<double> >& coords_sp) {
    if (coords_sp.size() != 3)
      throw std::invalid_argument("Error! Wrong dimension of field monitor points!");
    if (coords_sp[0].size() != coords_sp[1].size() || coords_sp[0].size() != coords_sp[2].size())
      throw std::invalid_argument("Error! Missing coordinates for field monitor points!");
    coords_sp_ = coords_sp;
    // for (int i = 0; i < coords_sp_[0].size(); ++i) {
    //   printf("%g, %g, %g\n", coords_sp_[0][i], coords_sp_[1][i], coords_sp_[2][i]);
    // }
  }  // end of void MultiLayerMie::SetFieldPointsSP(...)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetPEC(int layer_position) {
    isMieCalculated_ = false;
    if (layer_position < 0)
      throw std::invalid_argument("Error! Layers are numbered from 0!");
    PEC_layer_position_ = layer_position;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetMaxTermsNumber(int nmax) {    
    isMieCalculated_ = false;
    nmax_preset_ = nmax;
    //debug
    printf("Setting max terms: %d\n", nmax_preset_);
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::GenerateSizeParameter() {
    isMieCalculated_ = false;
    size_parameter_.clear();
    double radius = 0.0;
    for (auto width : target_width_) {
      radius += width;
      size_parameter_.push_back(2*PI_*radius/wavelength_);
    }
    for (auto width : coating_width_) {
      radius += width;
      size_parameter_.push_back(2*PI_*radius/wavelength_);
    }
    total_radius_ = radius;
  }  // end of void MultiLayerMie::GenerateSizeParameter();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::GenerateIndex() {
    isMieCalculated_ = false;
    index_.clear();
    for (auto index : target_index_) index_.push_back(index);
    for (auto index : coating_index_) index_.push_back(index);
  }  // end of void MultiLayerMie::GenerateIndex();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMie::GetTotalRadius() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    if (total_radius_ == 0) GenerateSizeParameter();
    return total_radius_;      
  }  // end of double MultiLayerMie::GetTotalRadius();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector< std::vector<double> >
  MultiLayerMie::GetSpectra(double from_WL, double to_WL, int samples) {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    std::vector< std::vector<double> > spectra;
    double step_WL = (to_WL - from_WL)/static_cast<double>(samples);
    double wavelength_backup = wavelength_;
    long fails = 0;
    for (double WL = from_WL; WL < to_WL; WL += step_WL) {
      wavelength_ = WL;
      try {
        RunMieCalculations();
      } catch(const std::invalid_argument& ia) {
        fails++;
        continue;
      }
      //printf("%3.1f ",WL);
      spectra.push_back(std::vector<double>({wavelength_, Qext_, Qsca_, Qabs_, Qbk_}));
    }  // end of for each WL in spectra
    printf("Spectrum has %li fails\n",fails);
    wavelength_ = wavelength_backup;
    return spectra;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::ClearTarget() {
    isMieCalculated_ = false;
    target_width_.clear();
    target_index_.clear();
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::ClearCoating() {
    isMieCalculated_ = false;
    coating_width_.clear();
    coating_index_.clear();
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::ClearLayers() {
    isMieCalculated_ = false;
    ClearTarget();
    ClearCoating();
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::ClearAllDesign() {
    isMieCalculated_ = false;
    ClearLayers();
    size_parameter_.clear();
    index_.clear();
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  //                         Computational core
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // Calculate Nstop - equation (17)
  //
  void MultiLayerMie::Nstop() {
    const double& xL = size_parameter_.back();
    if (xL <= 8) {
      nmax_ = round(xL + 4.0*pow(xL, 1.0/3.0) + 1);
    } else if (xL <= 4200) {
      nmax_ = round(xL + 4.05*pow(xL, 1.0/3.0) + 2);
    } else {
      nmax_ = round(xL + 4.0*pow(xL, 1.0/3.0) + 2);
    }    
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::Nmax(int first_layer) {
    int ri, riM1;
    const std::vector<double>& x = size_parameter_;
    const std::vector<std::complex<double> >& m = index_;
    Nstop();  // Set initial nmax_ value
    for (int i = first_layer; i < x.size(); i++) {
      if (i > PEC_layer_position_) 
        ri = round(std::abs(x[i]*m[i]));
      else 
        ri = 0;      
      nmax_ = std::max(nmax_, ri);
      // first layer is pec, if pec is present
      if ((i > first_layer) && ((i - 1) > PEC_layer_position_)) 
        riM1 = round(std::abs(x[i - 1]* m[i]));
      else 
        riM1 = 0;      
      nmax_ = std::max(nmax_, riM1);
    }
    nmax_ += 15;  // Final nmax_ value
  }
  //**********************************************************************************//
  // This function calculates the spherical Bessel (jn) and Hankel (h1n) functions    //
  // and their derivatives for a given complex value z. See pag. 87 B&H.              //
  //                                                                                  //
  // Input parameters:                                                                //
  //   z: Real argument to evaluate jn and h1n                                        //
  //   nmax_: Maximum number of terms to calculate jn and h1n                          //
  //                                                                                  //
  // Output parameters:                                                               //
  //   jn, h1n: Spherical Bessel and Hankel functions                                 //
  //   jnp, h1np: Derivatives of the spherical Bessel and Hankel functions            //
  //                                                                                  //
  // The implementation follows the algorithm by I.J. Thompson and A.R. Barnett,      //
  // Comp. Phys. Comm. 47 (1987) 245-257.                                             //
  //                                                                                  //
  // Complex spherical Bessel functions from n=0..nmax_ - 1 for z in the upper half      //
  // plane (Im(z) > -3).                                                              //
  //                                                                                  //
  //     j[n]   = j/n(z)                Regular solution: j[0]=sin(z)/z               //
  //     j'[n]  = d[j/n(z)]/dz                                                        //
  //     h1[n]  = h[0]/n(z)             Irregular Hankel function:                    //
  //     h1'[n] = d[h[0]/n(z)]/dz                h1[0] = j0(z) + i*y0(z)              //
  //                                                   = (sin(z)-i*cos(z))/z          //
  //                                                   = -i*exp(i*z)/z                //
  // Using complex CF1, and trigonometric forms for n=0 solutions.                    //
  //**********************************************************************************//
  void MultiLayerMie::sbesjh(std::complex<double> z,
                             std::vector<std::complex<double> >& jn,
                             std::vector<std::complex<double> >& jnp,
                             std::vector<std::complex<double> >& h1n,
                             std::vector<std::complex<double> >& h1np) {
    const int limit = 20000;
    const double accur = 1.0e-12;
    const double tm30 = 1e-30;

    double absc;
    std::complex<double> zi, w;
    std::complex<double> pl, f, b, d, c, del, jn0, jndb, h1nldb, h1nbdb;

    absc = std::abs(std::real(z)) + std::abs(std::imag(z));
    if ((absc < accur) || (std::imag(z) < -3.0)) {
      throw std::invalid_argument("TODO add error description for condition if ((absc < accur) || (std::imag(z) < -3.0))");
    }

    zi = 1.0/z;
    w = zi + zi;

    pl = double(nmax_)*zi;

    f = pl + zi;
    b = f + f + zi;
    d = 0.0;
    c = f;
    for (int n = 0; n < limit; n++) {
      d = b - d;
      c = b - 1.0/c;

      absc = std::abs(std::real(d)) + std::abs(std::imag(d));
      if (absc < tm30) {
        d = tm30;
      }

      absc = std::abs(std::real(c)) + std::abs(std::imag(c));
      if (absc < tm30) {
        c = tm30;
      }

      d = 1.0/d;
      del = d*c;
      f = f*del;
      b += w;

      absc = std::abs(std::real(del - 1.0)) + std::abs(std::imag(del - 1.0));

      if (absc < accur) {
        // We have obtained the desired accuracy
        break;
      }
    }

    if (absc > accur) {
      throw std::invalid_argument("We were not able to obtain the desired accuracy");
    }

    jn[nmax_ - 1] = tm30;
    jnp[nmax_ - 1] = f*jn[nmax_ - 1];

    // Downward recursion to n=0 (N.B.  Coulomb Functions)
    for (int n = nmax_ - 2; n >= 0; n--) {
      jn[n] = pl*jn[n + 1] + jnp[n + 1];
      jnp[n] = pl*jn[n] - jn[n + 1];
      pl = pl - zi;
    }

    // Calculate the n=0 Bessel Functions
    jn0 = zi*std::sin(z);
    h1n[0] = std::exp(std::complex<double>(0.0, 1.0)*z)*zi*(-std::complex<double>(0.0, 1.0));
    h1np[0] = h1n[0]*(std::complex<double>(0.0, 1.0) - zi);

    // Rescale j[n], j'[n], converting to spherical Bessel functions.
    // Recur   h1[n], h1'[n] as spherical Bessel functions.
    w = 1.0/jn[0];
    pl = zi;
    for (int n = 0; n < nmax_; n++) {
      jn[n] = jn0*(w*jn[n]);
      jnp[n] = jn0*(w*jnp[n]) - zi*jn[n];
      if (n != 0) {
        h1n[n] = (pl - zi)*h1n[n - 1] - h1np[n - 1];

        // check if hankel is increasing (upward stable)
        if (std::abs(h1n[n]) < std::abs(h1n[n - 1])) {
          jndb = z;
          h1nldb = h1n[n];
          h1nbdb = h1n[n - 1];
        }

        pl += zi;

        h1np[n] = -(pl*h1n[n]) + h1n[n - 1];
      }
    }
  }

  //**********************************************************************************//
  // This function calculates the spherical Bessel functions (bj and by) and the      //
  // logarithmic derivative (bd) for a given complex value z. See pag. 87 B&H.        //
  //                                                                                  //
  // Input parameters:                                                                //
  //   z: Complex argument to evaluate bj, by and bd                                  //
  //   nmax_: Maximum number of terms to calculate bj, by and bd                       //
  //                                                                                  //
  // Output parameters:                                                               //
  //   bj, by: Spherical Bessel functions                                             //
  //   bd: Logarithmic derivative                                                     //
  //**********************************************************************************//
  void MultiLayerMie::sphericalBessel(std::complex<double> z,
                                      std::vector<std::complex<double> >& bj,
                                      std::vector<std::complex<double> >& by,
                                      std::vector<std::complex<double> >& bd) {
    std::vector<std::complex<double> > jn(nmax_), jnp(nmax_), h1n(nmax_), h1np(nmax_);
    sbesjh(z, jn, jnp, h1n, h1np);

    for (int n = 0; n < nmax_; n++) {
      bj[n] = jn[n];
      by[n] = (h1n[n] - jn[n])/std::complex<double>(0.0, 1.0);
      bd[n] = jnp[n]/jn[n] + 1.0/z;
    }
    // std::complex<double> besselj_0 = std::sin(z)/z;
    // std::complex<double> bessely_0 = -std::cos(z)/z;
    // if (nmax_>0) {
    //   bj[0] = std::sin(z)/pow2(z)-std::cos(z)/z;  //bj1
    //   by[0] = std::cos(z)/pow2(z)-std::sin(z)/z;  //by1
    // }
    // if (nmax_>1) {
    //   bj[1] = bj[0]*3.0/z-besselj_0;//bj2
    //   by[1] = by[0]*3.0/z-bessely_0;//bj2
    // }
    // for (int n = 2; n < nmax_; n++) {
    //   bj[n] = (2.0*n - 1.0)/z*bj[n - 1] - bj[n];
    //   by[n] = (2.0*n - 1.0)/z*by[n - 1] - by[n];
    // }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // Calculate an - equation (5)
  std::complex<double> MultiLayerMie::calc_an(int n, double XL, std::complex<double> Ha, std::complex<double> mL,
                                              std::complex<double> PsiXL, std::complex<double> ZetaXL,
                                              std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1) {

    std::complex<double> Num = (Ha/mL + n/XL)*PsiXL - PsiXLM1;
    std::complex<double> Denom = (Ha/mL + n/XL)*ZetaXL - ZetaXLM1;

    return Num/Denom;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // Calculate bn - equation (6)
  std::complex<double> MultiLayerMie::calc_bn(int n, double XL, std::complex<double> Hb, std::complex<double> mL,
                                              std::complex<double> PsiXL, std::complex<double> ZetaXL,
                                              std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1) {

    std::complex<double> Num = (mL*Hb + n/XL)*PsiXL - PsiXLM1;
    std::complex<double> Denom = (mL*Hb + n/XL)*ZetaXL - ZetaXLM1;

    return Num/Denom;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // Calculates S1 - equation (25a)
  std::complex<double> MultiLayerMie::calc_S1(int n, std::complex<double> an, std::complex<double> bn,
                                              double Pi, double Tau) {
    return double(n + n + 1)*(Pi*an + Tau*bn)/double(n*n + n);
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // Calculates S2 - equation (25b) (it's the same as (25a), just switches Pi and Tau)
  std::complex<double> MultiLayerMie::calc_S2(int n, std::complex<double> an, std::complex<double> bn,
                                              double Pi, double Tau) {
    return calc_S1(n, an, bn, Tau, Pi);
  }
  //**********************************************************************************//
  // This function calculates the Riccati-Bessel functions (Psi and Zeta) for a       //
  // real argument (x).                                                               //
  // Equations (20a) - (21b)                                                          //
  //                                                                                  //
  // Input parameters:                                                                //
  //   x: Real argument to evaluate Psi and Zeta                                      //
  //   nmax: Maximum number of terms to calculate Psi and Zeta                        //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Psi, Zeta: Riccati-Bessel functions                                            //
  //**********************************************************************************//
  void MultiLayerMie::calcPsiZeta(std::complex<double> z,
                                  std::vector<std::complex<double> > D1,
                                  std::vector<std::complex<double> > D3,
                                  std::vector<std::complex<double> >& Psi,
                                  std::vector<std::complex<double> >& Zeta) {
    //Upward recurrence for Psi and Zeta - equations (20a) - (21b)
    //Psi[0] = std::complex<double>(std::sin(x), 0);
    std::complex<double> c_i(0.0, 1.0);
    Psi[0] = std::sin(z);
    //Zeta[0] = std::complex<double>(std::sin(x), -std::cos(x));
    Zeta[0] = std::sin(z) - c_i*std::cos(z);
    for (int n = 1; n <= nmax_; n++) {
      Psi[n] = Psi[n - 1]*(static_cast<double>(n)/z - D1[n - 1]);
      Zeta[n] = Zeta[n - 1]*(static_cast<double>(n)/z - D3[n - 1]);
    }

  }
  //**********************************************************************************//
  // Function CONFRA ported from MIEV0.f (Wiscombe,1979)
  // Ref. to NCAR Technical Notes, Wiscombe, 1979
  /*
c         Compute Bessel function ratio A-sub-N from its
c         continued fraction using Lentz method

c         ZINV = Reciprocal of argument of A


c    I N T E R N A L    V A R I A B L E S
c    ------------------------------------

c    CAK      Term in continued fraction expansion of A (Eq. R25)
c     a_k

c    CAPT     Factor used in Lentz iteration for A (Eq. R27)
c     T_k

c    CNUMER   Numerator   in capT  (Eq. R28A)
c     N_k
c    CDENOM   Denominator in capT  (Eq. R28B)
c     D_k

c    CDTD     Product of two successive denominators of capT factors
c                 (Eq. R34C)
c     xi_1

c    CNTN     Product of two successive numerators of capT factors
c                 (Eq. R34B)
c     xi_2

c    EPS1     Ill-conditioning criterion
c    EPS2     Convergence criterion

c    KK       Subscript k of cAk  (Eq. R25B)
c     k

c    KOUNT    Iteration counter (used to prevent infinite looping)

c    MAXIT    Max. allowed no. of iterations

c    MM + 1  and - 1, alternately
*/
  std::complex<double> MultiLayerMie::calcD1confra(const int N, const std::complex<double> z) {
  // NTMR -> nmax_ - 1  \\TODO nmax_ ?
    //int N = nmax_ - 1;
    int KK, KOUNT, MAXIT = 10000, MM;
    //    double EPS1=1.0e-2;
    double EPS2=1.0e-8;
    std::complex<double> CAK, CAPT, CDENOM, CDTD, CNTN, CNUMER;
    std::complex<double> one = std::complex<double>(1.0,0.0);
    std::complex<double> ZINV = one/z;
// c                                 ** Eq. R25a
    std::complex<double> CONFRA = static_cast<std::complex<double> >(N + 1)*ZINV;   //debug ZINV
    MM = - 1; 
    KK = 2*N +3; //debug 3
// c                                 ** Eq. R25b, k=2
    CAK    = static_cast<std::complex<double> >(MM*KK) * ZINV; //debug -3 ZINV
    CDENOM = CAK;
    CNUMER = CDENOM + one/CONFRA; //-3zinv+z
    KOUNT  = 1;
    //10 CONTINUE
    do {      ++KOUNT;
      if (KOUNT > MAXIT) {
        printf("re(%g):im(%g)\t\n", CONFRA.real(), CONFRA.imag());
        throw std::invalid_argument("ConFra--Iteration failed to converge!\n");
      }
      MM *= - 1;      KK += 2;  //debug  mm=1 kk=5
      CAK = static_cast<std::complex<double> >(MM*KK) * ZINV; //    ** Eq. R25b //debug 5zinv
     //  //c ** Eq. R32    Ill-conditioned case -- stride two terms instead of one
     //  if (std::abs(CNUMER/CAK) >= EPS1 ||  std::abs(CDENOM/CAK) >= EPS1) {
     //         //c                       ** Eq. R34
     //         CNTN   = CAK * CNUMER + 1.0;
     //         CDTD   = CAK * CDENOM + 1.0;
     //         CONFRA = (CNTN/CDTD) * CONFRA; // ** Eq. R33
     //         MM  *= - 1;        KK  += 2;
     //         CAK = static_cast<std::complex<double> >(MM*KK) * ZINV; // ** Eq. R25b
     //         //c                        ** Eq. R35
     //         CNUMER = CAK + CNUMER/CNTN;
     //         CDENOM = CAK + CDENOM/CDTD;
     //         ++KOUNT;
     //         //GO TO  10
     //         continue;
     // } else { //c                           *** Well-conditioned case
      {
        CAPT   = CNUMER/CDENOM; // ** Eq. R27 //debug (-3zinv + z)/(-3zinv)
        // printf("re(%g):im(%g)**\t", CAPT.real(), CAPT.imag());
       CONFRA = CAPT * CONFRA; // ** Eq. R26
       //if (N == 0) {output=true;printf(" re:");prn(CONFRA.real());printf(" im:"); prn(CONFRA.imag());output=false;};
       //c                                  ** Check for convergence; Eq. R31
       if (std::abs(CAPT.real() - 1.0) >= EPS2 ||  std::abs(CAPT.imag()) >= EPS2) {
//c                                        ** Eq. R30
         CNUMER = CAK + one/CNUMER;
         CDENOM = CAK + one/CDENOM;
         continue;
         //GO TO  10
       }  // end of if < eps2
      }
      break;
    } while(1);    
    //if (N == 0)  printf(" return confra for z=(%g,%g)\n", ZINV.real(), ZINV.imag());
    return CONFRA;
  }
  //**********************************************************************************//
  // This function calculates the logarithmic derivatives of the Riccati-Bessel       //
  // functions (D1 and D3) for a complex argument (z).                                //
  // Equations (16a), (16b) and (18a) - (18d)                                         //
  //                                                                                  //
  // Input parameters:                                                                //
  //   z: Complex argument to evaluate D1 and D3                                      //
  //   nmax_: Maximum number of terms to calculate D1 and D3                          //
  //                                                                                  //
  // Output parameters:                                                               //
  //   D1, D3: Logarithmic derivatives of the Riccati-Bessel functions                //
  //**********************************************************************************//
  void MultiLayerMie::calcD1D3(const std::complex<double> z,
                               std::vector<std::complex<double> >& D1,
                               std::vector<std::complex<double> >& D3) {
    // Downward recurrence for D1 - equations (16a) and (16b)
    D1[nmax_] = std::complex<double>(0.0, 0.0);
    //D1[nmax_] = calcD1confra(nmax_, z);
    const std::complex<double> zinv = std::complex<double>(1.0, 0.0)/z;
    
    // printf(" D:");prn((D1[nmax_]).real()); printf("\t diff:");
    // prn((D1[nmax_] + double(nmax_)*zinv).real());
    for (int n = nmax_; n > 0; n--) {
      D1[n - 1] = double(n)*zinv - 1.0/(D1[n] + double(n)*zinv);
      //D1[n - 1] = calcD1confra(n - 1, z);
      // printf(" D:");prn((D1[n - 1]).real()); printf("\t diff:");
      // prn((D1[n] + double(n)*zinv).real());
    }
    //     printf("\n\n"); iformat=0;
    if (std::abs(D1[0]) > 100000.0)
      throw std::invalid_argument
        ("Unstable D1! Please, try to change input parameters!\n");
    // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
    PsiZeta_[0] = 0.5*(1.0 - std::complex<double>(std::cos(2.0*z.real()), std::sin(2.0*z.real()))
                       *std::exp(-2.0*z.imag()));
    D3[0] = std::complex<double>(0.0, 1.0);
    for (int n = 1; n <= nmax_; n++) {
      PsiZeta_[n] = PsiZeta_[n - 1]*(static_cast<double>(n)*zinv - D1[n - 1])
        *(static_cast<double>(n)*zinv- D3[n - 1]);
      D3[n] = D1[n] + std::complex<double>(0.0, 1.0)/PsiZeta_[n];
    }
  }
  //**********************************************************************************//
  // This function calculates Pi and Tau for all values of Theta.                     //
  // Equations (26a) - (26c)                                                          //
  //                                                                                  //
  // Input parameters:                                                                //
  //   nmax_: Maximum number of terms to calculate Pi and Tau                         //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Pi, Tau: Angular functions Pi and Tau, as defined in equations (26a) - (26c)   //
  //**********************************************************************************//
  void MultiLayerMie::calcSinglePiTau(const double& costheta, std::vector<double>& Pi,
                                      std::vector<double>& Tau) {
    //****************************************************//
    // Equations (26a) - (26c)                            //
    //****************************************************//
    for (int n = 0; n < nmax_; n++) {
      if (n == 0) {
        // Initialize Pi and Tau
        Pi[n] = 1.0;
        Tau[n] = (n + 1)*costheta; 
      } else {
        // Calculate the actual values
        Pi[n] = ((n == 1) ? ((n + n + 1)*costheta*Pi[n - 1]/n)
                 : (((n + n + 1)*costheta*Pi[n - 1]
                     - (n + 1)*Pi[n - 2])/n));
        Tau[n] = (n + 1)*costheta*Pi[n] - (n + 2)*Pi[n - 1];
      }
    }
  }  // end of void MultiLayerMie::calcPiTau(...)
  void MultiLayerMie::calcAllPiTau(std::vector< std::vector<double> >& Pi,
                                   std::vector< std::vector<double> >& Tau) {
    std::vector<double> costheta(theta_.size(), 0.0);
    for (int t = 0; t < theta_.size(); t++) {
      costheta[t] = std::cos(theta_[t]);
    }
    // Do not join upper and lower 'for' to a single one!  It will slow
    // down the code!!! (For about 0.5-2.0% of runtime, it is probably
    // due to increased cache missing rate originated from the
    // recurrence in calcPiTau...)
    for (int t = 0; t < theta_.size(); t++) {
      calcSinglePiTau(costheta[t], Pi[t], Tau[t]);
      //calcSinglePiTau(std::cos(theta_[t]), Pi[t], Tau[t]); // It is slow!!
    }
  }  // end of void MultiLayerMie::calcAllPiTau(...)
  //**********************************************************************************//
  // This function calculates the scattering coefficients required to calculate       //
  // both the near- and far-field parameters.                                         //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it.             //
  //                                                                                  //
  // Output parameters:                                                               //
  //   an, bn: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  void MultiLayerMie::ScattCoeffs(std::vector<std::complex<double> >& an,
                                  std::vector<std::complex<double> >& bn) {
    const std::vector<double>& x = size_parameter_;
    const std::vector<std::complex<double> >& m = index_;
    const int& pl = PEC_layer_position_;
    const int L = index_.size();
    //************************************************************************//
    // Calculate the index of the first layer. It can be either 0 (default)   //
    // or the index of the outermost PEC layer. In the latter case all layers //
    // below the PEC are discarded.                                           //
    // ***********************************************************************//
    // TODO, is it possible for PEC to have a zero index? If yes than
    // is should be:
    // int fl = (pl > - 1) ? pl : 0;
    // This will give the same result, however, it corresponds the
    // logic - if there is PEC, than first layer is PEC.
    int fl = (pl > 0) ? pl : 0;
    if (nmax_ <= 0) Nmax(fl);

    std::complex<double> z1, z2;
    //**************************************************************************//
    // Note that since Fri, Nov 14, 2014 all arrays start from 0 (zero), which  //
    // means that index = layer number - 1 or index = n - 1. The only exception //
    // are the arrays for representing D1, D3 and Q because they need a value   //
    // for the index 0 (zero), hence it is important to consider this shift     //
    // between different arrays. The change was done to optimize memory usage.  //
    //**************************************************************************//
    // Allocate memory to the arrays
    std::vector<std::complex<double> > D1_mlxl(nmax_ + 1), D1_mlxlM1(nmax_ + 1),
      D3_mlxl(nmax_ + 1), D3_mlxlM1(nmax_ + 1);
    std::vector<std::vector<std::complex<double> > > Q(L), Ha(L), Hb(L);
    for (int l = 0; l < L; l++) {
      // D1_mlxl[l].resize(nmax_ + 1);
      // D1_mlxlM1[l].resize(nmax_ + 1);
      // D3_mlxl[l].resize(nmax_ + 1);
      // D3_mlxlM1[l].resize(nmax_ + 1);
      Q[l].resize(nmax_ + 1);
      Ha[l].resize(nmax_);
      Hb[l].resize(nmax_);
    }
    an.resize(nmax_);
    bn.resize(nmax_);
    PsiZeta_.resize(nmax_ + 1);
    std::vector<std::complex<double> > D1XL(nmax_ + 1), D3XL(nmax_ + 1), 
      PsiXL(nmax_ + 1), ZetaXL(nmax_ + 1);
    //*************************************************//
    // Calculate D1 and D3 for z1 in the first layer   //
    //*************************************************//
    if (fl == pl) {  // PEC layer
      for (int n = 0; n <= nmax_; n++) {
        D1_mlxl[n] = std::complex<double>(0.0, - 1.0);
        D3_mlxl[n] = std::complex<double>(0.0, 1.0);
      }
    } else { // Regular layer
      z1 = x[fl]* m[fl];
      // Calculate D1 and D3
      calcD1D3(z1, D1_mlxl, D3_mlxl);
    }
    // do { \
    //   ++iformat;\
    //   if (iformat%5 == 0) printf("%24.16e",z1.real());
    // } while (false);
    //******************************************************************//
    // Calculate Ha and Hb in the first layer - equations (7a) and (8a) //
    //******************************************************************//
    for (int n = 0; n < nmax_; n++) {
      Ha[fl][n] = D1_mlxl[n + 1];
      Hb[fl][n] = D1_mlxl[n + 1];
    }
    //*****************************************************//
    // Iteration from the second layer to the last one (L) //
    //*****************************************************//
    std::complex<double> Temp, Num, Denom;
    std::complex<double> G1, G2;
    for (int l = fl + 1; l < L; l++) {
      //************************************************************//
      //Calculate D1 and D3 for z1 and z2 in the layers fl + 1..L     //
      //************************************************************//
      z1 = x[l]*m[l];
      z2 = x[l - 1]*m[l];
      //Calculate D1 and D3 for z1
      calcD1D3(z1, D1_mlxl, D3_mlxl);
      //Calculate D1 and D3 for z2
      calcD1D3(z2, D1_mlxlM1, D3_mlxlM1);
      // prn(z1.real());
      // for (auto i : D1_mlxl) { prn(i.real());
      //   // prn(i.imag());
      //         } printf("\n");

      //*********************************************//
      //Calculate Q, Ha and Hb in the layers fl + 1..L //
      //*********************************************//
      // Upward recurrence for Q - equations (19a) and (19b)
      Num = std::exp(-2.0*(z1.imag() - z2.imag()))
        * std::complex<double>(std::cos(-2.0*z2.real()) - std::exp(-2.0*z2.imag()), std::sin(-2.0*z2.real()));
      Denom = std::complex<double>(std::cos(-2.0*z1.real()) - std::exp(-2.0*z1.imag()), std::sin(-2.0*z1.real()));
      Q[l][0] = Num/Denom;
      for (int n = 1; n <= nmax_; n++) {
        Num = (z1*D1_mlxl[n] + double(n))*(double(n) - z1*D3_mlxl[n - 1]);
        Denom = (z2*D1_mlxlM1[n] + double(n))*(double(n) - z2*D3_mlxlM1[n - 1]);
        Q[l][n] = ((pow2(x[l - 1]/x[l])* Q[l][n - 1])*Num)/Denom;
      }
      // Upward recurrence for Ha and Hb - equations (7b), (8b) and (12) - (15)
      for (int n = 1; n <= nmax_; n++) {
        //Ha
        if ((l - 1) == pl) { // The layer below the current one is a PEC layer
          G1 = -D1_mlxlM1[n];
          G2 = -D3_mlxlM1[n];
        } else {
          G1 = (m[l]*Ha[l - 1][n - 1]) - (m[l - 1]*D1_mlxlM1[n]);
          G2 = (m[l]*Ha[l - 1][n - 1]) - (m[l - 1]*D3_mlxlM1[n]);
        }  // end of if PEC
        Temp = Q[l][n]*G1;
        Num = (G2*D1_mlxl[n]) - (Temp*D3_mlxl[n]);
        Denom = G2 - Temp;
        Ha[l][n - 1] = Num/Denom;
        //Hb
        if ((l - 1) == pl) { // The layer below the current one is a PEC layer
          G1 = Hb[l - 1][n - 1];
          G2 = Hb[l - 1][n - 1];
        } else {
          G1 = (m[l - 1]*Hb[l - 1][n - 1]) - (m[l]*D1_mlxlM1[n]);
          G2 = (m[l - 1]*Hb[l - 1][n - 1]) - (m[l]*D3_mlxlM1[n]);
        }  // end of if PEC

        Temp = Q[l][n]*G1;
        Num = (G2*D1_mlxl[n]) - (Temp* D3_mlxl[n]);
        Denom = (G2- Temp);
        Hb[l][n - 1] = (Num/ Denom);
      }  // end of for Ha and Hb terms
    }  // end of for layers iteration
    //**************************************//
    //Calculate D1, D3, Psi and Zeta for XL //
    //**************************************//
    // Calculate D1XL and D3XL
    calcD1D3(x[L - 1],  D1XL, D3XL);
    //printf("%5.20f\n",Ha[L - 1][0].real());
    // Calculate PsiXL and ZetaXL
    calcPsiZeta(x[L - 1], D1XL, D3XL, PsiXL, ZetaXL);
    //*********************************************************************//
    // Finally, we calculate the scattering coefficients (an and bn) and   //
    // the angular functions (Pi and Tau). Note that for these arrays the  //
    // first layer is 0 (zero), in future versions all arrays will follow  //
    // this convention to save memory. (13 Nov, 2014)                      //
    //*********************************************************************//
    for (int n = 0; n < nmax_; n++) {
      //********************************************************************//
      //Expressions for calculating an and bn coefficients are not valid if //
      //there is only one PEC layer (ie, for a simple PEC sphere).          //
      //********************************************************************//
      if (pl < (L - 1)) {
        an[n] = calc_an(n + 1, x[L - 1], Ha[L - 1][n], m[L - 1], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
        bn[n] = calc_bn(n + 1, x[L - 1], Hb[L - 1][n], m[L - 1], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
      } else {
        an[n] = calc_an(n + 1, x[L - 1], std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
        bn[n] = PsiXL[n + 1]/ZetaXL[n + 1];
      }
    }  // end of for an and bn terms
  }  // end of void MultiLayerMie::ScattCoeffs(...)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::InitMieCalculations() {
    isMieCalculated_ = false;
    // Initialize the scattering parameters
    Qext_ = 0;
    Qsca_ = 0;
    Qabs_ = 0;
    Qbk_ = 0;
    Qpr_ = 0;
    asymmetry_factor_ = 0;
    albedo_ = 0;
    Qsca_ch_.clear();
    Qext_ch_.clear();
    Qabs_ch_.clear();
    Qbk_ch_.clear();
    Qpr_ch_.clear();
    Qsca_ch_.resize(nmax_ - 1);
    Qext_ch_.resize(nmax_ - 1);
    Qabs_ch_.resize(nmax_ - 1);
    Qbk_ch_.resize(nmax_ - 1);
    Qpr_ch_.resize(nmax_ - 1);
    Qsca_ch_norm_.resize(nmax_ - 1);
    Qext_ch_norm_.resize(nmax_ - 1);
    Qabs_ch_norm_.resize(nmax_ - 1);
    Qbk_ch_norm_.resize(nmax_ - 1);
    Qpr_ch_norm_.resize(nmax_ - 1);
    // Initialize the scattering amplitudes
    std::vector<std::complex<double> >        tmp1(theta_.size(),std::complex<double>(0.0, 0.0));
    S1_.swap(tmp1);
    S2_ = S1_;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::ConvertToSP() {
    isMieCalculated_ = false;
    if (target_width_.size() + coating_width_.size() == 0)
      return;  // Nothing to convert, we suppose that SP was set directly
    GenerateSizeParameter();
    GenerateIndex();
    if (size_parameter_.size() != index_.size())
      throw std::invalid_argument("Ivalid conversion of width to size parameter units!/n");
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  //**********************************************************************************//
  // This function calculates the actual scattering parameters and amplitudes         //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send - 1                          //
  //   x: Array containing the size parameters of the layers [0..L - 1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L - 1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //   nmax_: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to - 1 and the function will calculate it              //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  void MultiLayerMie::RunMieCalculations() {
    isMieCalculated_ = false;
    ConvertToSP();
    nmax_ = nmax_preset_;
    if (size_parameter_.size() != index_.size())
      throw std::invalid_argument("Each size parameter should have only one index!");
    if (size_parameter_.size() == 0)
      throw std::invalid_argument("Initialize model first!");
    const std::vector<double>& x = size_parameter_;
    // Calculate scattering coefficients
    ScattCoeffs(an_, bn_);

    // std::vector< std::vector<double> > Pi(nmax_), Tau(nmax_);
    std::vector< std::vector<double> > Pi, Tau;
    Pi.resize(theta_.size());
    Tau.resize(theta_.size());
    for (int i =0; i< theta_.size(); ++i) {
      Pi[i].resize(nmax_);
      Tau[i].resize(nmax_);
    }
    calcAllPiTau(Pi, Tau);    
    InitMieCalculations(); //
    std::complex<double> Qbktmp(0.0, 0.0);
    std::vector< std::complex<double> > Qbktmp_ch(nmax_ - 1, Qbktmp);
    // By using downward recurrence we avoid loss of precision due to float rounding errors
    // See: https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
    //      http://en.wikipedia.org/wiki/Loss_of_significance
    for (int i = nmax_ - 2; i >= 0; i--) {
      const int n = i + 1;
      // Equation (27)
      Qext_ch_norm_[i] = (an_[i].real() + bn_[i].real());
      Qext_ch_[i] = (n + n + 1.0)*Qext_ch_norm_[i];
      //Qext_ch_[i] = (n + n + 1)*(an_[i].real() + bn_[i].real());
      Qext_ += Qext_ch_[i];
      // Equation (28)
      Qsca_ch_norm_[i] = (an_[i].real()*an_[i].real() + an_[i].imag()*an_[i].imag()
                          + bn_[i].real()*bn_[i].real() + bn_[i].imag()*bn_[i].imag());
      Qsca_ch_[i] = (n + n + 1.0)*Qsca_ch_norm_[i];
      Qsca_ += Qsca_ch_[i];
      // Qsca_ch_[i] += (n + n + 1)*(an_[i].real()*an_[i].real() + an_[i].imag()*an_[i].imag()
      //                             + bn_[i].real()*bn_[i].real() + bn_[i].imag()*bn_[i].imag());

      // Equation (29) TODO We must check carefully this equation. If we
      // remove the typecast to double then the result changes. Which is
      // the correct one??? Ovidio (2014/12/10) With cast ratio will
      // give double, without cast (n + n + 1)/(n*(n + 1)) will be
      // rounded to integer. Tig (2015/02/24)
      Qpr_ch_[i]=((n*(n + 2)/(n + 1))*((an_[i]*std::conj(an_[n]) + bn_[i]*std::conj(bn_[n])).real())
               + ((double)(n + n + 1)/(n*(n + 1)))*(an_[i]*std::conj(bn_[i])).real());
      Qpr_ += Qpr_ch_[i];
      // Equation (33)      
      Qbktmp_ch[i] = (double)(n + n + 1)*(1 - 2*(n % 2))*(an_[i]- bn_[i]);
      Qbktmp += Qbktmp_ch[i];
      // Calculate the scattering amplitudes (S1 and S2)    //
      // Equations (25a) - (25b)                            //
      for (int t = 0; t < theta_.size(); t++) {
        S1_[t] += calc_S1(n, an_[i], bn_[i], Pi[t][i], Tau[t][i]);
        S2_[t] += calc_S2(n, an_[i], bn_[i], Pi[t][i], Tau[t][i]);
      }
    }
    double x2 = pow2(x.back());
    Qext_ = 2.0*(Qext_)/x2;                                 // Equation (27)
    for (double& Q : Qext_ch_) Q = 2.0*Q/x2;
    Qsca_ = 2.0*(Qsca_)/x2;                                 // Equation (28)
    for (double& Q : Qsca_ch_) Q = 2.0*Q/x2;
    //for (double& Q : Qsca_ch_norm_) Q = 2.0*Q/x2;
    Qpr_ = Qext_ - 4.0*(Qpr_)/x2;                           // Equation (29)
    for (int i = 0; i < nmax_ - 1; ++i) Qpr_ch_[i] = Qext_ch_[i] - 4.0*Qpr_ch_[i]/x2;

    Qabs_ = Qext_ - Qsca_;                                // Equation (30)
    for (int i = 0; i < nmax_ - 1; ++i) {
      Qabs_ch_[i] = Qext_ch_[i] - Qsca_ch_[i];
      Qabs_ch_norm_[i] = Qext_ch_norm_[i] - Qsca_ch_norm_[i];
    }
    
    albedo_ = Qsca_/Qext_;                              // Equation (31)
    asymmetry_factor_ = (Qext_ - Qpr_)/Qsca_;                          // Equation (32)

    Qbk_ = (Qbktmp.real()*Qbktmp.real() + Qbktmp.imag()*Qbktmp.imag())/x2;    // Equation (33)

    isMieCalculated_ = true;
    nmax_used_ = nmax_;
    printf("Run Mie result: Qext = %g, Qsca = %g, Qabs = %g, Qbk = %g \n",
               GetQext(), GetQsca(), GetQabs(), GetQbk());
    //return nmax;
  }
  

  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::ScattCoeffsLayerdInit() {
    const int L = index_.size();
    // we need to fill
    // std::vector< std::vector<std::complex<double> > > al_n_, bl_n_, cl_n_, dl_n_;
    //     for n = [0..nmax_) and for l=[L..0)
    // TODO: to decrease cache miss outer loop is with n and inner with reversed l
    // at the moment outer is forward l and inner in n
    al_n_.resize(L + 1);
    bl_n_.resize(L + 1);
    cl_n_.resize(L + 1);
    dl_n_.resize(L + 1);
    for (auto& element:al_n_) element.resize(nmax_);
    for (auto& element:bl_n_) element.resize(nmax_);
    for (auto& element:cl_n_) element.resize(nmax_);
    for (auto& element:dl_n_) element.resize(nmax_);
    std::complex<double> c_one(1.0, 0.0);
    std::complex<double> c_zero(0.0, 0.0);
    // Yang, paragraph under eq. A3
    // a^(L + 1)_n = a_n, d^(L + 1) = 1 ...
    for (int i = 0; i < nmax_; ++i) {
      al_n_[L][i] = an_[i];
      bl_n_[L][i] = bn_[i];
      cl_n_[L][i] = c_one;
      dl_n_[L][i] = c_one;
      if (i<3) printf(" (%g) ", std::abs(an_[i]));
    }

  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::ScattCoeffsLayerd() {
    if (!isMieCalculated_)
      throw std::invalid_argument("(ScattCoeffsLayerd) You should run calculations first!");
    ScattCoeffsLayerdInit();
    const int L = index_.size();
    std::vector<std::complex<double> > z(L), z1(L);
    for (int i = 0; i < L - 1; ++i) {
      z[i]  =size_parameter_[i]*index_[i];
      z1[i]=size_parameter_[i]*index_[i + 1];
    }
    z[L - 1]  =size_parameter_[L - 1]*index_[L - 1];
    z1[L - 1]  =size_parameter_[L - 1];
    std::vector< std::vector<std::complex<double> > > D1z(L), D1z1(L), D3z(L), D3z1(L);
    std::vector< std::vector<std::complex<double> > > Psiz(L), Psiz1(L), Zetaz(L), Zetaz1(L);
    for (int l = 0; l < L; ++l) {
      D1z[l].resize(nmax_ + 1);
      D1z1[l].resize(nmax_ + 1);
      D3z[l].resize(nmax_ + 1);
      D3z1[l].resize(nmax_ + 1);
      Psiz[l].resize(nmax_ + 1);
      Psiz1[l].resize(nmax_ + 1);
      Zetaz[l].resize(nmax_ + 1);
      Zetaz1[l].resize(nmax_ + 1);
    }
    for (int l = 0; l < L; ++l) {
      calcD1D3(z[l],D1z[l],D3z[l]);
      calcD1D3(z1[l],D1z1[l],D3z1[l]);
      calcPsiZeta(z[l],D1z[l],D3z[l], Psiz[l],Zetaz[l]);
      calcPsiZeta(z1[l],D1z1[l],D3z1[l], Psiz1[l],Zetaz1[l]);
    }
    auto& m = index_;
    std::vector< std::complex<double> > m1(L);
    for (int l = 0; l < L - 1; ++l) m1[l] = m[l + 1];
    m1[L - 1] = std::complex<double> (1.0, 0.0);
    // for (auto zz : m) printf ("m[i]=%g \n\n ", zz.real());
    for (int l = L - 1; l >= 0; --l) {
      for (int n = 0; n < nmax_; ++n) {
        // al_n
        auto denom = m1[l]*Zetaz[l][n + 1] * (D1z[l][n + 1] - D3z[l][n + 1]);
        al_n_[l][n] = D1z[l][n + 1]* m1[l]
          *(al_n_[l + 1][n]*Zetaz1[l][n + 1] - dl_n_[l + 1][n]*Psiz1[l][n + 1])
          - m[l]*(-D1z1[l][n + 1]*dl_n_[l + 1][n]*Psiz1[l][n + 1]
                  +D3z1[l][n + 1]*al_n_[l + 1][n]*Zetaz1[l][n + 1]);
        al_n_[l][n] /= denom;
        // if (n<2) printf("denom[%d][%d]:%g \n", l, n,
        //                   std::abs(Psiz[l][n + 1]));
        // dl_n
        denom = m1[l]*Psiz[l][n + 1] * (D1z[l][n + 1] - D3z[l][n + 1]);
        dl_n_[l][n] = D3z[l][n + 1]*m1[l]
          *(al_n_[l + 1][n]*Zetaz1[l][n + 1] - dl_n_[l + 1][n]*Psiz1[l][n + 1])
          - m[l]*(-D1z1[l][n + 1]*dl_n_[l + 1][n]*Psiz1[l][n + 1]
                  +D3z1[l][n + 1]*al_n_[l + 1][n]*Zetaz1[l][n + 1]);
        dl_n_[l][n] /= denom;
        // bl_n
        denom = m1[l]*Zetaz[l][n + 1] * (D1z[l][n + 1] - D3z[l][n + 1]);
        bl_n_[l][n] = D1z[l][n + 1]* m[l]
          *(bl_n_[l + 1][n]*Zetaz1[l][n + 1] - cl_n_[l + 1][n]*Psiz1[l][n + 1])
          - m1[l]*(-D1z1[l][n + 1]*cl_n_[l + 1][n]*Psiz1[l][n + 1]
                  +D3z1[l][n + 1]*bl_n_[l + 1][n]*Zetaz1[l][n + 1]);
        bl_n_[l][n] /= denom;
        // cl_n
        denom = m1[l]*Psiz[l][n + 1] * (D1z[l][n + 1] - D3z[l][n + 1]);
        cl_n_[l][n] = D3z[l][n + 1]*m[l]
          *(bl_n_[l + 1][n]*Zetaz1[l][n + 1] - cl_n_[l + 1][n]*Psiz1[l][n + 1])
          - m1[l]*(-D1z1[l][n + 1]*cl_n_[l + 1][n]*Psiz1[l][n + 1]
                  +D3z1[l][n + 1]*bl_n_[l + 1][n]*Zetaz1[l][n + 1]);
        cl_n_[l][n] /= denom;   
      }  // end of all n
    }  // end of for all l
    // Check the result and change  an__0 and bn__0 for exact zero
    for (int n = 0; n < nmax_; ++n) {
      if (std::abs(al_n_[0][n]) < 1e-10) al_n_[0][n] = 0.0;
      else throw std::invalid_argument("Unstable calculation of a__0_n!");
      if (std::abs(bl_n_[0][n]) < 1e-10) bl_n_[0][n] = 0.0;
      else throw std::invalid_argument("Unstable calculation of b__0_n!");
    }

    // for (int l = 0; l < L; ++l) {
    //   printf("l=%d --> ", l);
    //   for (int n = 0; n < nmax_ + 1; ++n) {
    //         if (n < 20) continue;
    //         printf("n=%d --> D1zn=%g, D3zn=%g, D1zn=%g, D3zn=%g || ",
    //                n,
    //                D1z[l][n].real(), D3z[l][n].real(),
    //                D1z1[l][n].real(), D3z1[l][n].real());
    //   }
    //   printf("\n\n");
    // }
    // for (int l = 0; l < L; ++l) {
    //   printf("l=%d --> ", l);
    //   for (int n = 0; n < nmax_ + 1; ++n) {
    //         printf("n=%d --> D1zn=%g, D3zn=%g, D1zn=%g, D3zn=%g || ",
    //                n,
    //                D1z[l][n].real(), D3z[l][n].real(),
    //                D1z1[l][n].real(), D3z1[l][n].real());
    //   }
    //   printf("\n\n");
    // }
    for (int i = 0; i < L + 1; ++i) {
      printf("Layer =%d ---> ", i);
      for (int n = 0; n < nmax_; ++n) {
            //        if (n < 20) continue;
            printf(" || n=%d --> a=%g,%g b=%g,%g c=%g,%g d=%g,%g",
                   n,
                   al_n_[i][n].real(), al_n_[i][n].imag(),
                   bl_n_[i][n].real(), bl_n_[i][n].imag(),
                   cl_n_[i][n].real(), cl_n_[i][n].imag(),
                   dl_n_[i][n].real(), dl_n_[i][n].imag());
      }
      printf("\n\n");
    }
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  // external scattering field = incident + scattered
  // BH p.92 (4.37), 94 (4.45), 95 (4.50)
  // assume: medium is non-absorbing; refim = 0; Uabs = 0

  void MultiLayerMie::fieldExt(const double Rho, const double Phi, const double Theta, const  std::vector<double>& Pi, const std::vector<double>& Tau, std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H)  {
    
    std::complex<double> c_zero(0.0, 0.0), c_i(0.0, 1.0);
    std::vector<std::complex<double> > vm3o1n(3), vm3e1n(3), vn3o1n(3), vn3e1n(3);
    std::vector<std::complex<double> > Ei(3,c_zero), Hi(3,c_zero), Es(3,c_zero), Hs(3,c_zero);
    std::vector<std::complex<double> > bj(nmax_ + 1), by(nmax_ + 1), bd(nmax_ + 1);
    // Calculate spherical Bessel and Hankel functions
    printf("##########  layer OUT ############\n");
    sphericalBessel(Rho,bj, by, bd);    
    for (int n = 0; n < nmax_; n++) {
      double rn = static_cast<double>(n + 1);
      std::complex<double> zn = bj[n + 1] + c_i*by[n + 1];
      // using BH 4.12 and 4.50
      std::complex<double> xxip = Rho*(bj[n] + c_i*by[n]) - rn*zn;
      
      using std::sin;
      using std::cos;
      vm3o1n[0] = c_zero;
      vm3o1n[1] = cos(Phi)*Pi[n]*zn;
      vm3o1n[2] = -sin(Phi)*Tau[n]*zn;
      vm3e1n[0] = c_zero;
      vm3e1n[1] = -sin(Phi)*Pi[n]*zn;
      vm3e1n[2] = -cos(Phi)*Tau[n]*zn;
      vn3o1n[0] = sin(Phi)*rn*(rn + 1.0)*sin(Theta)*Pi[n]*zn/Rho;
      vn3o1n[1] = sin(Phi)*Tau[n]*xxip/Rho;
      vn3o1n[2] = cos(Phi)*Pi[n]*xxip/Rho;
      vn3e1n[0] = cos(Phi)*rn*(rn + 1.0)*sin(Theta)*Pi[n]*zn/Rho;
      vn3e1n[1] = cos(Phi)*Tau[n]*xxip/Rho;
      vn3e1n[2] = -sin(Phi)*Pi[n]*xxip/Rho;
      
      // scattered field: BH p.94 (4.45)
      std::complex<double> encap = std::pow(c_i, rn)*(2.0*rn + 1.0)/(rn*rn + rn);
      for (int i = 0; i < 3; i++) {
        Es[i] = Es[i] + encap*(c_i*an_[n]*vn3e1n[i] - bn_[n]*vm3o1n[i]);
        Hs[i] = Hs[i] + encap*(c_i*bn_[n]*vn3o1n[i] + an_[n]*vm3e1n[i]);
        //if (n<3) printf(" E[%d]=%g ", i,std::abs(Es[i]));
        if (n<3) printf(" !!=%d=== %g ", i,std::abs(Es[i]));
      }
    }
    
    // incident E field: BH p.89 (4.21); cf. p.92 (4.37), p.93 (4.38)
    // basis unit vectors = er, etheta, ephi
    std::complex<double> eifac = std::exp(std::complex<double>(0.0, Rho*std::cos(Theta)));
    {
      using std::sin;
      using std::cos;
      Ei[0] = eifac*sin(Theta)*cos(Phi);
      Ei[1] = eifac*cos(Theta)*cos(Phi);
      Ei[2] = -eifac*sin(Phi);
    }

    // magnetic field
    double hffact = 1.0/(cc_*mu_);
    for (int i = 0; i < 3; i++) {
      Hs[i] = hffact*Hs[i];
    }
    
    // incident H field: BH p.26 (2.43), p.89 (4.21)
    std::complex<double> hffacta = hffact;
    std::complex<double> hifac = eifac*hffacta;
    {
      using std::sin;
      using std::cos;
      Hi[0] = hifac*sin(Theta)*sin(Phi);
      Hi[1] = hifac*cos(Theta)*sin(Phi);
      Hi[2] = hifac*cos(Phi);
    }
    
    for (int i = 0; i < 3; i++) {
      // electric field E [V m - 1] = EF*E0
      E[i] = Ei[i] + Es[i];
      H[i] = Hi[i] + Hs[i];
      // printf("ext E[%d]=%g",i,std::abs(E[i]));
    }
   }  // end of void fieldExt(...)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::fieldInt(const double Rho, const double Phi, const double Theta, const  std::vector<double>& Pi, const std::vector<double>& Tau, std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H)  {
    // printf("field int Qext = %g, Qsca = %g, Qabs = %g, Qbk = %g, \n",
    //            GetQext(), GetQsca(), GetQabs(), GetQbk());
    
    std::complex<double> c_zero(0.0, 0.0), c_i(0.0, 1.0), c_one(1.0, 0.0);
    std::vector<std::complex<double> > vm3o1n(3), vm3e1n(3), vn3o1n(3), vn3e1n(3);
    std::vector<std::complex<double> > vm1o1n(3), vm1e1n(3), vn1o1n(3), vn1e1n(3);
    std::vector<std::complex<double> > El(3,c_zero),Ei(3,c_zero), Hl(3,c_zero);
    std::vector<std::complex<double> > bj(nmax_ + 1), by(nmax_ + 1), bd(nmax_ + 1);
    int layer=0;  // layer number
    std::complex<double> index_l;
    for (int i = 0; i < size_parameter_.size() - 1; ++i) {
      if (size_parameter_[i] < Rho && Rho <= size_parameter_[i + 1]) {
        layer=i;
      }
    }
    if (Rho > size_parameter_.back()) {
      layer = size_parameter_.size();
      index_l = c_one; 
    } else {
      index_l = index_[layer]; 
    }
   
    std::complex<double> bessel_arg = Rho*index_l;
    std::complex<double>& rh = bessel_arg;
    std::complex<double> besselj_1 = std::sin(rh)/pow2(rh)-std::cos(rh)/rh;
    printf("bessel arg = %g,%g   index=%g,%g   besselj[1]=%g,%g\n", bessel_arg.real(), bessel_arg.imag(), index_l.real(), index_l.imag(), besselj_1.real(), besselj_1.imag());
    const int& l = layer;
    printf("##########  layer %d ############\n",l);
    // Calculate spherical Bessel and Hankel functions
    sphericalBessel(bessel_arg,bj, by, bd);    
    printf("besselj[1]=%g,%g\n", bj[1].real(), bj[1].imag());
    printf("bessely[1]=%g,%g\n", by[1].real(), by[1].imag());
    for (int n = 0; n < nmax_; n++) {
      double rn = static_cast<double>(n + 1);
      std::complex<double> znm1 = bj[n] + c_i*by[n];
      std::complex<double> zn = bj[n + 1] + c_i*by[n + 1];
      //if (n<3) printf("\nbesselh = %g,%g", zn.real(), zn.imag()); //!
      // using BH 4.12 and 4.50
      std::complex<double> xxip = Rho*(bj[n] + c_i*by[n]) - rn*zn;
      //if (n<3) printf("\nxxip = %g,%g", xxip.real(), xxip.imag()); //!
      
      using std::sin;
      using std::cos;
      vm3o1n[0] = c_zero;
      vm3o1n[1] = cos(Phi)*Pi[n]*zn;
      vm3o1n[2] = -sin(Phi)*Tau[n]*zn;
      // if (n<3)  printf("\nRE  vm3o1n[0]%g   vm3o1n[1]%g    vm3o1n[2]%g   \nIM vm3o1n[0]%g   vm3o1n[1]%g    vm3o1n[2]%g",
      //              vm3o1n[0].real(),   vm3o1n[1].real(),    vm3o1n[2].real(),
      //              vm3o1n[0].imag(),   vm3o1n[1].imag(),    vm3o1n[2].imag());
      vm3e1n[0] = c_zero;
      vm3e1n[1] = -sin(Phi)*Pi[n]*zn;
      vm3e1n[2] = -cos(Phi)*Tau[n]*zn;
      vn3o1n[0] = sin(Phi)*rn*(rn + 1.0)*sin(Theta)*Pi[n]*zn/Rho;
      vn3o1n[1] = sin(Phi)*Tau[n]*xxip/Rho;
      vn3o1n[2] = cos(Phi)*Pi[n]*xxip/Rho;
      vn3e1n[0] = cos(Phi)*rn*(rn + 1.0)*sin(Theta)*Pi[n]*zn/Rho;
      vn3e1n[1] = cos(Phi)*Tau[n]*xxip/Rho;
      vn3e1n[2] = -sin(Phi)*Pi[n]*xxip/Rho;
      // if (n<3)  printf("\nRE  vn3e1n[0]%g   vn3e1n[1]%g    vn3e1n[2]%g   \nIM vn3e1n[0]%g   vn3e1n[1]%g    vn3e1n[2]%g",
      //              vn3e1n[0].real(),   vn3e1n[1].real(),    vn3e1n[2].real(),
      //              vn3e1n[0].imag(),   vn3e1n[1].imag(),    vn3e1n[2].imag());
      
      znm1 = bj[n];
      zn = bj[n + 1];
      // znm1 = (bj[n] + c_i*by[n]).real();
      // zn = (bj[n + 1] + c_i*by[n + 1]).real();
      xxip = Rho*(bj[n]) - rn*zn;
      if (n<3)printf("\nbesselj = %g,%g", zn.real(), zn.imag()); //!
      vm1o1n[0] = c_zero;
      vm1o1n[1] = cos(Phi)*Pi[n]*zn;
      vm1o1n[2] = -sin(Phi)*Tau[n]*zn;
      vm1e1n[0] = c_zero;
      vm1e1n[1] = -sin(Phi)*Pi[n]*zn;
      vm1e1n[2] = -cos(Phi)*Tau[n]*zn;
      vn1o1n[0] = sin(Phi)*rn*(rn + 1.0)*sin(Theta)*Pi[n]*zn/Rho;
      vn1o1n[1] = sin(Phi)*Tau[n]*xxip/Rho;
      vn1o1n[2] = cos(Phi)*Pi[n]*xxip/Rho;
      // if (n<3) printf("\nvn1o1n[2](%g) = cos(Phi)(%g)*Pi[n](%g)*xxip(%g)/Rho(%g)",
      //                       std::abs(vn1o1n[2]), cos(Phi),Pi[n],std::abs(xxip),Rho);
      vn1e1n[0] = cos(Phi)*rn*(rn + 1.0)*sin(Theta)*Pi[n]*zn/Rho;
      vn1e1n[1] = cos(Phi)*Tau[n]*xxip/Rho;
      vn1e1n[2] = -sin(Phi)*Pi[n]*xxip/Rho;
      // if (n<3)  printf("\nRE  vm3o1n[0]%g   vm3o1n[1]%g    vm3o1n[2]%g   \nIM vm3o1n[0]%g   vm3o1n[1]%g    vm3o1n[2]%g",
      //              vm3o1n[0].real(),   vm3o1n[1].real(),    vm3o1n[2].real(),
      //              vm3o1n[0].imag(),   vm3o1n[1].imag(),    vm3o1n[2].imag());
      
      // scattered field: BH p.94 (4.45)
      std::complex<double> encap = std::pow(c_i, rn)*(2.0*rn + 1.0)/(rn*rn + rn);
      // if (n<3) printf("\n===== n=%d ======\n",n);
      for (int i = 0; i < 3; i++) {
        // if (n<3 && i==0) printf("\nn=%d",n);
        // if (n<3) printf("\nbefore !El[%d]=%g,%g! ", i, El[i].real(), El[i].imag());
        Ei[i] = encap*(
                       cl_n_[l][n]*vm1o1n[i] - c_i*dl_n_[l][n]*vn1e1n[i]
          + c_i*al_n_[l][n]*vn3e1n[i] - bl_n_[l][n]*vm3o1n[i]
);
        El[i] = El[i] + encap*(cl_n_[l][n]*vm1o1n[i] - c_i*dl_n_[l][n]*vn1e1n[i]
                               + c_i*al_n_[l][n]*vn3e1n[i] - bl_n_[l][n]*vm3o1n[i]);
        Hl[i] = Hl[i] + encap*(-dl_n_[l][n]*vm1e1n[i] - c_i*cl_n_[l][n]*vn1o1n[i]
                               + c_i*bl_n_[l][n]*vn3o1n[i] + al_n_[l][n]*vm3e1n[i]);
        // printf("\n !Ei[%d]=%g,%g! ", i, Ei[i].real(), Ei[i].imag());
        // if (n<3) printf("\n !El[%d]=%g,%g! ", i, El[i].real(), El[i].imag());
        // //printf(" ===%d=== %g ", i,std::abs(cl_n_[l][n]*vm1o1n[i] - c_i*dl_n_[l][n]*vn1e1n[i]));
        // if (n<3) printf(" ===%d=== %g ", i,std::abs(//-dl_n_[l][n]*vm1e1n[i] 
        //                                             //- c_i*cl_n_[l][n]*
        //                                             vn1o1n[i]
        //                                             // + c_i*bl_n_[l][n]*vn3o1n[i]
        //                                             // + al_n_[l][n]*vm3e1n[i]
        //                      ));
        // if (n<3) printf(" --- Ei[%d]=%g! ", i,std::abs(encap*(vm1o1n[i] - c_i*vn1e1n[i])));

      }
      //if (n<3) printf(" bj=%g \n", std::abs(bj[n]));
    }  // end of for all n
    
    // magnetic field
    double hffact = 1.0/(cc_*mu_);
    for (int i = 0; i < 3; i++) {
      Hl[i] = hffact*Hl[i];
    }
    
    for (int i = 0; i < 3; i++) {
      // electric field E [V m - 1] = EF*E0
      E[i] = El[i];
      H[i] = Hl[i];
      printf("\n !El[%d]=%g,%g! ", i, El[i].real(), El[i].imag());
      //printf(" E[%d]=%g",i,std::abs(El[i]));
    }
   }  // end of void fieldExt(...)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //

  //**********************************************************************************//
  // This function calculates complex electric and magnetic field in the surroundings //
  // and inside (TODO) the particle.                                                  //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send 0 (zero)                    //
  //   x: Array containing the size parameters of the layers [0..L - 1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L - 1]     //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to 0 (zero) and the function will calculate it.       //
  //   ncoord: Number of coordinate points                                            //
  //   Coords: Array containing all coordinates where the complex electric and        //
  //           magnetic fields will be calculated                                     //
  //                                                                                  //
  // Output parameters:                                                               //
  //   E, H: Complex electric and magnetic field at the provided coordinates          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  void MultiLayerMie::RunFieldCalculations() {
    // Calculate scattering coefficients an_ and bn_
    RunMieCalculations();
    //nmax_=10;
    ScattCoeffsLayerd();

    std::vector<double> Pi(nmax_), Tau(nmax_);
    long total_points = coords_sp_[0].size();
    E_field_.resize(total_points);
    H_field_.resize(total_points);
    for (auto& f:E_field_) f.resize(3);
    for (auto& f:H_field_) f.resize(3);

    for (int point = 0; point < total_points; ++point) {
      const double& Xp = coords_sp_[0][point];
      const double& Yp = coords_sp_[1][point];
      const double& Zp = coords_sp_[2][point];
      printf("X=%g, Y=%g, Z=%g\n", Xp, Yp, Zp);
      // Convert to spherical coordinates
      double Rho, Phi, Theta;
      Rho = std::sqrt(pow2(Xp) + pow2(Yp) + pow2(Zp));
      // printf("Rho=%g\n", Rho);
      // Avoid convergence problems due to Rho too small
      if (Rho < 1e-10) Rho = 1e-10;
      // If Rho=0 then Theta is undefined. Just set it to zero to avoid problems
      if (Rho == 0.0) Theta = 0.0;
      else Theta = std::acos(Zp/Rho);
      // printf("Theta=%g\n", Theta);
      // If Xp=Yp=0 then Phi is undefined. Just set it to zero to zero to avoid problems
      if (Xp == 0.0 && Yp == 0.0) Phi = 0.0;
      else Phi = std::acos(Xp/std::sqrt(pow2(Xp) + pow2(Yp)));
      // printf("Phi=%g\n", Phi);

      calcSinglePiTau(std::cos(Theta), Pi, Tau);     
      //*******************************************************//
      // external scattering field = incident + scattered      //
      // BH p.92 (4.37), 94 (4.45), 95 (4.50)                  //
      // assume: medium is non-absorbing; refim = 0; Uabs = 0  //
      //*******************************************************//
      // This array contains the fields in spherical coordinates
      std::vector<std::complex<double> > Es(3), Hs(3);
      const double outer_size = size_parameter_.back();
      // Firstly the easiest case: the field outside the particle
      printf("rho=%g, outer=%g  ", Rho, outer_size);
      if (Rho >= outer_size) {
        fieldExt(Rho, Phi, Theta, Pi, Tau, Es, Hs);
        printf("\nFin E ext: %g,%g,%g   Rho=%g\n", std::abs(Es[0]), std::abs(Es[1]),std::abs(Es[2]), Rho);
      } else {
        fieldInt(Rho, Phi, Theta, Pi, Tau, Es, Hs);      
        printf("\nFin E int: %g,%g,%g   Rho=%g\n", std::abs(Es[0]), std::abs(Es[1]),std::abs(Es[2]), Rho);
      }
      std::complex<double>& Ex = E_field_[point][0];
      std::complex<double>& Ey = E_field_[point][1];
      std::complex<double>& Ez = E_field_[point][2];
      std::complex<double>& Hx = H_field_[point][0];
      std::complex<double>& Hy = H_field_[point][1];
      std::complex<double>& Hz = H_field_[point][2];
      //Now, convert the fields back to cartesian coordinates
      {
        using std::sin;
        using std::cos;
        Ex = sin(Theta)*cos(Phi)*Es[0] + cos(Theta)*cos(Phi)*Es[1] - sin(Phi)*Es[2];
        Ey = sin(Theta)*sin(Phi)*Es[0] + cos(Theta)*sin(Phi)*Es[1] + cos(Phi)*Es[2];
        Ez = cos(Theta)*Es[0] - sin(Theta)*Es[1];
      
        Hx = sin(Theta)*cos(Phi)*Hs[0] + cos(Theta)*cos(Phi)*Hs[1] - sin(Phi)*Hs[2];
        Hy = sin(Theta)*sin(Phi)*Hs[0] + cos(Theta)*sin(Phi)*Hs[1] + cos(Phi)*Hs[2];
        Hz = cos(Theta)*Hs[0] - sin(Theta)*Hs[1];
      }
      printf("Cart E: %g,%g,%g   Rho=%g\n", std::abs(Ex), std::abs(Ey),std::abs(Ez),
             Rho);
    }  // end of for all field coordinates
    
  }  //  end of   void MultiLayerMie::RunFieldCalculations()

}  // end of namespace nmie
