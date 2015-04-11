//**********************************************************************************//
//    Copyright (C) 2009-2015  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2015  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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

//**********************************************************************************//
// This class implements the algorithm for a multilayered sphere described by:      //
//    [1] W. Yang, "Improved recursive algorithm for light scattering by a          //
//        multilayered sphere,” Applied Optics, vol. 42, Mar. 2003, pp. 1710-1720.  //
//                                                                                  //
// You can find the description of all the used equations in:                       //
//    [2] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//                                                                                  //
// Hereinafter all equations numbers refer to [2]                                   //
//**********************************************************************************//
#include "nmie.h"
#include <array>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace nmie {
  //helpers
  template<class T> inline T pow2(const T value) {return value*value;}

  int round(double x) {
    return x >= 0 ? (int)(x + 0.5):(int)(x - 0.5);
  }


  //**********************************************************************************//
  // This function emulates a C call to calculate the actual scattering parameters    //
  // and amplitudes.                                                                  //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it              //
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
  int nMie(const int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {

    if (x.size() != L || m.size() != L)
        throw std::invalid_argument("Declared number of layers do not fit x and m!");
    if (Theta.size() != nTheta)
        throw std::invalid_argument("Declared number of sample for Theta is not correct!");
    try {
      MultiLayerMie multi_layer_mie;
      multi_layer_mie.SetLayersSize(x);
      multi_layer_mie.SetLayersIndex(m);
      multi_layer_mie.SetAngles(Theta);

      multi_layer_mie.RunMieCalculation();

      *Qext = multi_layer_mie.GetQext();
      *Qsca = multi_layer_mie.GetQsca();
      *Qabs = multi_layer_mie.GetQabs();
      *Qbk = multi_layer_mie.GetQbk();
      *Qpr = multi_layer_mie.GetQpr();
      *g = multi_layer_mie.GetAsymmetryFactor();
      *Albedo = multi_layer_mie.GetAlbedo();
      S1 = multi_layer_mie.GetS1();
      S2 = multi_layer_mie.GetS2();
    } catch(const std::invalid_argument& ia) {
      // Will catch if  multi_layer_mie fails or other errors.
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      throw std::invalid_argument(ia);
      return -1;
    }

    return 0;
  }

  //**********************************************************************************//
  // This function is just a wrapper to call the full 'nMie' function with fewer      //
  // parameters, it is here mainly for compatibility with older versions of the       //
  // program. Also, you can use it if you neither have a PEC layer nor want to define //
  // any limit for the maximum number of terms.                                       //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
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
  int nMie(const int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
    return nmie::nMie(L, -1, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }


  //**********************************************************************************//
  // This function is just a wrapper to call the full 'nMie' function with fewer      //
  // parameters, it is useful if you want to include a PEC layer but not a limit      //
  // for the maximum number of terms.                                                 //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
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
  int nMie(const int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
    return nmie::nMie(L, pl, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }

  //**********************************************************************************//
  // This function is just a wrapper to call the full 'nMie' function with fewer      //
  // parameters, it is useful if you want to include a limit for the maximum number   //
  // of terms but not a PEC layer.                                                    //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it              //
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
  int nMie(const int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
    return nmie::nMie(L, -1, x, m, nTheta, Theta, nmax, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }


  //**********************************************************************************//
  // This function emulates a C call to calculate complex electric and magnetic field //
  // in the surroundings and inside (TODO) the particle.                              //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send 0 (zero)                    //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
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
  int nField(const int L, const int pl, const std::vector<double>& x, const std::vector<std::complex<double> >& m, const int nmax, const int ncoord, const std::vector<double>& Xp_vec, const std::vector<double>& Yp_vec, const std::vector<double>& Zp_vec, std::vector<std::vector<std::complex<double> > >& E, std::vector<std::vector<std::complex<double> > >& H) {
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
      //multi_layer_mie.SetPECLayer(pl);
      multi_layer_mie.SetLayersSize(x);
      multi_layer_mie.SetLayersIndex(m);
      multi_layer_mie.SetFieldCoords({Xp_vec, Yp_vec, Zp_vec});
      multi_layer_mie.RunFieldCalculation();
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
  // Returns previously calculated Qext                                     //
  // ********************************************************************** //
  double MultiLayerMie::GetQext() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qext_;
  }


  // ********************************************************************** //
  // Returns previously calculated Qabs                                     //
  // ********************************************************************** //
  double MultiLayerMie::GetQabs() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qabs_;
  }


  // ********************************************************************** //
  // Returns previously calculated Qsca                                     //
  // ********************************************************************** //
  double MultiLayerMie::GetQsca() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qsca_;
  }


  // ********************************************************************** //
  // Returns previously calculated Qbk                                      //
  // ********************************************************************** //
  double MultiLayerMie::GetQbk() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qbk_;
  }


  // ********************************************************************** //
  // Returns previously calculated Qpr                                      //
  // ********************************************************************** //
  double MultiLayerMie::GetQpr() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return Qpr_;
  }


  // ********************************************************************** //
  // Returns previously calculated assymetry factor                         //
  // ********************************************************************** //
  double MultiLayerMie::GetAsymmetryFactor() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return asymmetry_factor_;
  }


  // ********************************************************************** //
  // Returns previously calculated Albedo                                   //
  // ********************************************************************** //
  double MultiLayerMie::GetAlbedo() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return albedo_;
  }


  // ********************************************************************** //
  // Returns previously calculated S1                                       //
  // ********************************************************************** //
  std::vector<std::complex<double> > MultiLayerMie::GetS1() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return S1_;
  }


  // ********************************************************************** //
  // Returns previously calculated S2                                       //
  // ********************************************************************** //
  std::vector<std::complex<double> > MultiLayerMie::GetS2() {
    if (!isMieCalculated_)
      throw std::invalid_argument("You should run calculations before result request!");
    return S2_;
  }


  // ********************************************************************** //
  // Modify scattering (theta) angles                                       //
  // ********************************************************************** //
  void MultiLayerMie::SetAngles(const std::vector<double>& angles) {
    areIntCoeffsCalc_ = false;
    areExtCoeffsCalc_ = false;
    isMieCalculated_ = false;
    theta_ = angles;
  }


  // ********************************************************************** //
  // Modify size of all layers                                             //
  // ********************************************************************** //
  void MultiLayerMie::SetLayersSize(const std::vector<double>& layer_size) {
    areIntCoeffsCalc_ = false;
    areExtCoeffsCalc_ = false;
    isMieCalculated_ = false;
    size_param_.clear();
    double prev_layer_size = 0.0;
    for (auto curr_layer_size : layer_size) {
      if (curr_layer_size <= 0.0)
        throw std::invalid_argument("Size parameter should be positive!");
      if (prev_layer_size > curr_layer_size)
        throw std::invalid_argument
          ("Size parameter for next layer should be larger than the previous one!");
      prev_layer_size = curr_layer_size;
      size_param_.push_back(curr_layer_size);
    }
  }


  // ********************************************************************** //
  // Modify refractive index of all layers                                  //
  // ********************************************************************** //
  void MultiLayerMie::SetLayersIndex(const std::vector< std::complex<double> >& index) {
    areIntCoeffsCalc_ = false;
    areExtCoeffsCalc_ = false;
    isMieCalculated_ = false;
    refr_index_ = index;
  }


  // ********************************************************************** //
  // Modify coordinates for field calculation                               //
  // ********************************************************************** //
  void MultiLayerMie::SetFieldCoords(const std::vector< std::vector<double> >& coords) {
    if (coords.size() != 3)
      throw std::invalid_argument("Error! Wrong dimension of field monitor points!");
    if (coords[0].size() != coords[1].size() || coords[0].size() != coords[2].size())
      throw std::invalid_argument("Error! Missing coordinates for field monitor points!");
    coords_ = coords;
  }


  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetPECLayer(int layer_position) {
    areIntCoeffsCalc_ = false;
    areExtCoeffsCalc_ = false;
    isMieCalculated_ = false;
    if (layer_position < 0)
      throw std::invalid_argument("Error! Layers are numbered from 0!");
    PEC_layer_position_ = layer_position;
  }


  // ********************************************************************** //
  // Set maximun number of terms to be used                                 //
  // ********************************************************************** //
  void MultiLayerMie::SetMaxTerms(int nmax) {
    areIntCoeffsCalc_ = false;
    areExtCoeffsCalc_ = false;
    isMieCalculated_ = false;
    nmax_preset_ = nmax;
  }


  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMie::GetSizeParameter() {
    if (size_param_.size() > 0)
      return size_param_.back();
    else
      return 0;
  }


  // ********************************************************************** //
  // Clear layer information                                                //
  // ********************************************************************** //
  void MultiLayerMie::ClearLayers() {
    areIntCoeffsCalc_ = false;
    areExtCoeffsCalc_ = false;
    isMieCalculated_ = false;
    size_param_.clear();
    refr_index_.clear();
  }


  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  //                         Computational core
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //


  // ********************************************************************** //
  // Calculate calcNstop - equation (17)                                        //
  // ********************************************************************** //
  void MultiLayerMie::calcNstop() {
    const double& xL = size_param_.back();
    if (xL <= 8) {
      nmax_ = round(xL + 4.0*pow(xL, 1.0/3.0) + 1);
    } else if (xL <= 4200) {
      nmax_ = round(xL + 4.05*pow(xL, 1.0/3.0) + 2);
    } else {
      nmax_ = round(xL + 4.0*pow(xL, 1.0/3.0) + 2);
    }
  }


  // ********************************************************************** //
  // Maximum number of terms required for the calculation                   //
  // ********************************************************************** //
  void MultiLayerMie::calcNmax(int first_layer) {
    int ri, riM1;
    const std::vector<double>& x = size_param_;
    const std::vector<std::complex<double> >& m = refr_index_;
    calcNstop();  // Set initial nmax_ value
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
  //   nmax_: Maximum number of terms to calculate jn and h1n                         //
  //                                                                                  //
  // Output parameters:                                                               //
  //   jn, h1n: Spherical Bessel and Hankel functions                                 //
  //   jnp, h1np: Derivatives of the spherical Bessel and Hankel functions            //
  //                                                                                  //
  // The implementation follows the algorithm by I.J. Thompson and A.R. Barnett,      //
  // Comp. Phys. Comm. 47 (1987) 245-257.                                             //
  //                                                                                  //
  // Complex spherical Bessel functions from n=0..nmax_-1 for z in the upper half     //
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
  }


  // ********************************************************************** //
  // Calculate an - equation (5)                                            //
  // ********************************************************************** //
  std::complex<double> MultiLayerMie::calc_an(int n, double XL, std::complex<double> Ha, std::complex<double> mL,
                                              std::complex<double> PsiXL, std::complex<double> ZetaXL,
                                              std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1) {

    std::complex<double> Num = (Ha/mL + n/XL)*PsiXL - PsiXLM1;
    std::complex<double> Denom = (Ha/mL + n/XL)*ZetaXL - ZetaXLM1;

    return Num/Denom;
  }


  // ********************************************************************** //
  // Calculate bn - equation (6)                                            //
  // ********************************************************************** //
  std::complex<double> MultiLayerMie::calc_bn(int n, double XL, std::complex<double> Hb, std::complex<double> mL,
                                              std::complex<double> PsiXL, std::complex<double> ZetaXL,
                                              std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1) {

    std::complex<double> Num = (mL*Hb + n/XL)*PsiXL - PsiXLM1;
    std::complex<double> Denom = (mL*Hb + n/XL)*ZetaXL - ZetaXLM1;

    return Num/Denom;
  }


  // ********************************************************************** //
  // Calculates S1 - equation (25a)                                         //
  // ********************************************************************** //
  std::complex<double> MultiLayerMie::calc_S1(int n, std::complex<double> an, std::complex<double> bn,
                                              double Pi, double Tau) {
    return double(n + n + 1)*(Pi*an + Tau*bn)/double(n*n + n);
  }


  // ********************************************************************** //
  // Calculates S2 - equation (25b) (it's the same as (25a), just switches  //
  // Pi and Tau)                                                            //
  // ********************************************************************** //
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
    std::complex<double> c_i(0.0, 1.0);
    Psi[0] = std::sin(z);
    Zeta[0] = std::sin(z) - c_i*std::cos(z);
    for (int n = 1; n <= nmax_; n++) {
      Psi[n] = Psi[n - 1]*(static_cast<double>(n)/z - D1[n - 1]);
      Zeta[n] = Zeta[n - 1]*(static_cast<double>(n)/z - D3[n - 1]);
    }

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
    const std::complex<double> zinv = std::complex<double>(1.0, 0.0)/z;

    for (int n = nmax_; n > 0; n--) {
      D1[n - 1] = double(n)*zinv - 1.0/(D1[n] + double(n)*zinv);
    }

    if (std::abs(D1[0]) > 100000.0)
      throw std::invalid_argument("Unstable D1! Please, try to change input parameters!\n");

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
  // This function calculates Pi and Tau for a given value of cos(Theta).             //
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
  void MultiLayerMie::calcPiTau(const double& costheta,
                                std::vector<double>& Pi, std::vector<double>& Tau) {

    int n;
    //****************************************************//
    // Equations (26a) - (26c)                            //
    //****************************************************//
    // Initialize Pi and Tau
    Pi[0] = 1.0;
    Tau[0] = costheta;
    // Calculate the actual values
    if (nmax_ > 1) {
      Pi[1] = 3*costheta*Pi[0];
      Tau[1] = 2*costheta*Pi[1] - 3*Pi[0];
      for (n = 2; n < nmax_; n++) {
        Pi[n] = ((n + n + 1)*costheta*Pi[n - 1] - (n + 1)*Pi[n - 2])/n;
        Tau[n] = (n + 1)*costheta*Pi[n] - (n + 2)*Pi[n - 1];
      }
    }
  }  // end of MultiLayerMie::calcPiTau(...)

  void MultiLayerMie::calcSpherHarm(const double Rho, const double Phi, const double Theta,
                                    const std::complex<double>& zn, const std::complex<double>& deriv,
                                    const double& Pi, const double& Tau, const double& rn,
                                    std::vector<std::complex<double> >& Mo1n, std::vector<std::complex<double> >& Me1n, 
                                    std::vector<std::complex<double> >& No1n, std::vector<std::complex<double> >& Ne1n) {

    // using BH 4.12 and 4.50
    std::complex<double> c_zero(0.0, 0.0);

    using std::sin;
    using std::cos;
    Mo1n[0] = c_zero;
    Mo1n[1] = cos(Phi)*Pi*zn;
    Mo1n[2] = -sin(Phi)*Tau*zn;
    Me1n[0] = c_zero;
    Me1n[1] = -sin(Phi)*Pi*zn;
    Me1n[2] = -cos(Phi)*Tau*zn;
    No1n[0] = sin(Phi)*rn*(rn + 1.0)*sin(Theta)*Pi*zn/Rho;
    No1n[1] = sin(Phi)*Tau*deriv/Rho;
    No1n[2] = cos(Phi)*Pi*deriv/Rho;
    Ne1n[0] = cos(Phi)*rn*(rn + 1.0)*sin(Theta)*Pi*zn/Rho;
    Ne1n[1] = cos(Phi)*Tau*deriv/Rho;
    Ne1n[2] = -sin(Phi)*Pi*deriv/Rho;
  }  // end of MultiLayerMie::calcSpherHarm(...)


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
  void MultiLayerMie::ExtScattCoeffs() {

    areExtCoeffsCalc_ = false;

    const std::vector<double>& x = size_param_;
    const std::vector<std::complex<double> >& m = refr_index_;
    const int& pl = PEC_layer_position_;
    const int L = refr_index_.size();

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
    // Well, I followed the logic: First layer is always zero unless it has
    // an upper PEC layer.
    int fl = (pl > 0) ? pl : 0;
    if (nmax_preset_ <= 0) calcNmax(fl);
    else nmax_ = nmax_preset_;

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
      Q[l].resize(nmax_ + 1);
      Ha[l].resize(nmax_);
      Hb[l].resize(nmax_);
    }

    an_.resize(nmax_);
    bn_.resize(nmax_);
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

      //*********************************************//
      //Calculate Q, Ha and Hb in the layers fl + 1..L //
      //*********************************************//
      // Upward recurrence for Q - equations (19a) and (19b)
      Num = std::exp(-2.0*(z1.imag() - z2.imag()))
           *std::complex<double>(std::cos(-2.0*z2.real()) - std::exp(-2.0*z2.imag()), std::sin(-2.0*z2.real()));
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
    calcD1D3(x[L - 1], D1XL, D3XL);

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
        an_[n] = calc_an(n + 1, x[L - 1], Ha[L - 1][n], m[L - 1], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
        bn_[n] = calc_bn(n + 1, x[L - 1], Hb[L - 1][n], m[L - 1], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
      } else {
        an_[n] = calc_an(n + 1, x[L - 1], std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
        bn_[n] = PsiXL[n + 1]/ZetaXL[n + 1];
      }
    }  // end of for an and bn terms
    areExtCoeffsCalc_ = true;
  }  // end of void MultiLayerMie::ExtScattCoeffs(...)

  //**********************************************************************************//
  // This function calculates the actual scattering parameters and amplitudes         //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //   nmax_: Maximum number of multipolar expansion terms to be used for the         //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it              //
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
  void MultiLayerMie::RunMieCalculation() {
    if (size_param_.size() != refr_index_.size())
      throw std::invalid_argument("Each size parameter should have only one index!");
    if (size_param_.size() == 0)
      throw std::invalid_argument("Initialize model first!");

    const std::vector<double>& x = size_param_;

    areIntCoeffsCalc_ = false;
    areExtCoeffsCalc_ = false;
    isMieCalculated_ = false;

    // Calculate scattering coefficients
    ExtScattCoeffs();

//    for (int i = 0; i < nmax_; i++) {
//      printf("a[%i] = %g, %g; b[%i] = %g, %g\n", i, an_[i].real(), an_[i].imag(), i, bn_[i].real(), bn_[i].imag());
//    }

    if (!areExtCoeffsCalc_)
      throw std::invalid_argument("Calculation of scattering coefficients failed!");

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
    std::vector<std::complex<double> > tmp1(theta_.size(),std::complex<double>(0.0, 0.0));
    S1_.swap(tmp1);
    S2_ = S1_;

    std::vector<double> Pi(nmax_), Tau(nmax_);

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
        calcPiTau(std::cos(theta_[t]), Pi, Tau);

        S1_[t] += calc_S1(n, an_[i], bn_[i], Pi[i], Tau[i]);
        S2_[t] += calc_S2(n, an_[i], bn_[i], Pi[i], Tau[i]);
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
  }


  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::IntScattCoeffs() {
    if (!areExtCoeffsCalc_)
      throw std::invalid_argument("(IntScattCoeffs) You should calculate external coefficients first!");

    areIntCoeffsCalc_ = false;

    std::complex<double> c_one(1.0, 0.0);
    std::complex<double> c_zero(0.0, 0.0);

    const int L = refr_index_.size();

    // we need to fill
    // std::vector< std::vector<std::complex<double> > > anl_, bnl_, cnl_, dnl_;
    //     for n = [0..nmax_) and for l=[L..0)
    // TODO: to decrease cache miss outer loop is with n and inner with reversed l
    // at the moment outer is forward l and inner in n
    anl_.resize(L + 1);
    bnl_.resize(L + 1);
    cnl_.resize(L + 1);
    dnl_.resize(L + 1);
    for (auto& element:anl_) element.resize(nmax_);
    for (auto& element:bnl_) element.resize(nmax_);
    for (auto& element:cnl_) element.resize(nmax_);
    for (auto& element:dnl_) element.resize(nmax_);

    // Yang, paragraph under eq. A3
    // a^(L + 1)_n = a_n, d^(L + 1) = 1 ...
    for (int i = 0; i < nmax_; ++i) {
      anl_[L][i] = an_[i];
      bnl_[L][i] = bn_[i];
      cnl_[L][i] = c_one;
      dnl_[L][i] = c_one;
    }

    std::vector<std::complex<double> > z(L), z1(L);
    for (int i = 0; i < L - 1; ++i) {
      z[i] = size_param_[i]*refr_index_[i];
      z1[i] = size_param_[i]*refr_index_[i + 1];
    }
    z[L - 1] = size_param_[L - 1]*refr_index_[L - 1];
    z1[L - 1] = size_param_[L - 1];
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
      calcD1D3(z[l], D1z[l], D3z[l]);
      calcD1D3(z1[l], D1z1[l], D3z1[l]);
      calcPsiZeta(z[l], D1z[l], D3z[l], Psiz[l], Zetaz[l]);
      calcPsiZeta(z1[l], D1z1[l], D3z1[l], Psiz1[l], Zetaz1[l]);
    }
    auto& m = refr_index_;
    std::vector< std::complex<double> > m1(L);

    for (int l = 0; l < L - 1; ++l) m1[l] = m[l + 1];
    m1[L - 1] = std::complex<double> (1.0, 0.0);

    // for (auto zz : m) printf ("m[i]=%g \n\n ", zz.real());
    for (int l = L - 1; l >= 0; l--) {
      for (int n = nmax_ - 2; n >= 0; n--) {
        auto denomZeta = m1[l]*Zetaz[l][n + 1]*(D1z[l][n + 1] - D3z[l][n + 1]);
        auto denomPsi = m1[l]*Psiz[l][n + 1]*(D1z[l][n + 1] - D3z[l][n + 1]);

        auto T1 = anl_[l + 1][n]*Zetaz1[l][n + 1] - dnl_[l + 1][n]*Psiz1[l][n + 1];
        auto T2 = bnl_[l + 1][n]*Zetaz1[l][n + 1] - cnl_[l + 1][n]*Psiz1[l][n + 1];

        auto T3 = -D1z1[l][n + 1]*dnl_[l + 1][n]*Psiz1[l][n + 1] + D3z1[l][n + 1]*anl_[l + 1][n]*Zetaz1[l][n + 1];
        auto T4 = -D1z1[l][n + 1]*cnl_[l + 1][n]*Psiz1[l][n + 1] + D3z1[l][n + 1]*bnl_[l + 1][n]*Zetaz1[l][n + 1];

        // anl
        anl_[l][n] = (D1z[l][n + 1]*m1[l]*T1 - m[l]*T3)/denomZeta;
        // bnl
        bnl_[l][n] = (D1z[l][n + 1]*m[l]*T2 - m1[l]*T4)/denomZeta;
        // cnl
        cnl_[l][n] = (D3z[l][n + 1]*m[l]*T2 - m1[l]*T4)/denomPsi;
        // dnl
        dnl_[l][n] = (D3z[l][n + 1]*m1[l]*T1 - m[l]*T3)/denomPsi;
      }  // end of all n
    }  // end of all l

    // Check the result and change  an__0 and bn__0 for exact zero
    for (int n = 0; n < nmax_; ++n) {
      if (std::abs(anl_[0][n]) < 1e-10) anl_[0][n] = 0.0;
      else throw std::invalid_argument("Unstable calculation of a__0_n!");
      if (std::abs(bnl_[0][n]) < 1e-10) bnl_[0][n] = 0.0;
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
    //for (int i = 0; i < L + 1; ++i) {
    //  printf("Layer =%d ---> \n", i);
    //  for (int n = 0; n < nmax_; ++n) {
    //                if (n < 20) continue;
    //        printf(" || n=%d --> a=%g,%g b=%g,%g c=%g,%g d=%g,%g\n",
    //               n,
    //               anl_[i][n].real(), anl_[i][n].imag(),
    //               bnl_[i][n].real(), bnl_[i][n].imag(),
    //               cnl_[i][n].real(), cnl_[i][n].imag(),
    //               dnl_[i][n].real(), dnl_[i][n].imag());
    //  }
    //  printf("\n\n");
    //}
    areIntCoeffsCalc_ = true;
  }
  // ********************************************************************** //
  // ********************************************************************** //


  // ********************************************************************** //
  // external scattering field = incident + scattered                       //
  // BH p.92 (4.37), 94 (4.45), 95 (4.50)                                   //
  // assume: medium is non-absorbing; refim = 0; Uabs = 0                   //
  // ********************************************************************** //

  void MultiLayerMie::fieldExt(const double Rho, const double Phi, const double Theta, const  std::vector<double>& Pi, const std::vector<double>& Tau, std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H)  {

    std::complex<double> c_zero(0.0, 0.0), c_i(0.0, 1.0);
    std::vector<std::complex<double> > M3o1n(3), M3e1n(3), N3o1n(3), N3e1n(3);
    std::vector<std::complex<double> > Ei(3, c_zero), Hi(3, c_zero), Es(3, c_zero), Hs(3, c_zero);
    std::vector<std::complex<double> > bj(nmax_ + 1), by(nmax_ + 1), bd(nmax_ + 1);

    // Calculate spherical Bessel and Hankel functions
    sphericalBessel(Rho, bj, by, bd);

    //printf("##########  layer OUT ############\n");
    for (int n = 0; n < nmax_; n++) {
      int n1 = n + 1;
      double rn = static_cast<double>(n1);

      std::complex<double> zn = bj[n1] + c_i*by[n1];
      // using BH 4.12 and 4.50
      std::complex<double> deriv = Rho*(bj[n] + c_i*by[n]) - rn*zn;

      calcSpherHarm(Rho, Phi, Theta, zn, deriv, Pi[n], Tau[n], rn, M3o1n, M3e1n, N3o1n, N3e1n);

      // scattered field: BH p.94 (4.45)
      std::complex<double> En = std::pow(c_i, rn)*(2.0*rn + 1.0)/(rn*rn + rn);
      for (int i = 0; i < 3; i++) {
        Es[i] = Es[i] + En*(c_i*an_[n]*N3e1n[i] - bn_[n]*M3o1n[i]);
        Hs[i] = Hs[i] + En*(c_i*bn_[n]*N3o1n[i] + an_[n]*M3e1n[i]);
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
   }  // end of MultiLayerMie::fieldExt(...)


  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::fieldInt(const double Rho, const double Phi, const double Theta, const  std::vector<double>& Pi, const std::vector<double>& Tau, std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H)  {

    std::complex<double> c_zero(0.0, 0.0), c_i(0.0, 1.0), c_one(1.0, 0.0);
    std::vector<std::complex<double> > M3o1n(3), M3e1n(3), N3o1n(3), N3e1n(3);
    std::vector<std::complex<double> > M1o1n(3), M1e1n(3), N1o1n(3), N1e1n(3);
    std::vector<std::complex<double> > El(3, c_zero), Ei(3, c_zero), Hl(3, c_zero);
    std::vector<std::complex<double> > bj(nmax_ + 1), by(nmax_ + 1), bd(nmax_ + 1);

    int l = 0;  // Layer number
    std::complex<double> ml;

    if (Rho > size_param_.back()) {
      l = size_param_.size();
      ml = c_one;
    } else {
      for (int i = 0; i < size_param_.size() - 1; ++i) {
        if (size_param_[i] < Rho && Rho <= size_param_[i + 1]) {
          l = i;
        }
      }
      ml = refr_index_[l];
    }

    // Calculate spherical Bessel and Hankel functions
    sphericalBessel(Rho*ml, bj, by, bd);

    for (int n = nmax_ - 1; n >= 0; n--) {
      int n1 = n + 1;
      double rn = static_cast<double>(n1);

      // using BH 4.12 and 4.50

      std::complex<double> zn = bj[n1];
      std::complex<double> deriv = Rho*bj[n] - rn*zn;

      calcSpherHarm(Rho, Phi, Theta, zn, deriv, Pi[n], Tau[n], rn, M1o1n, M1e1n, N1o1n, N1e1n);

      zn = bj[n1] + c_i*by[n1];
      deriv = Rho*(bj[n] + c_i*by[n]) - rn*zn;

      calcSpherHarm(Rho, Phi, Theta, zn, deriv, Pi[n], Tau[n], rn, M3o1n, M3e1n, N3o1n, N3e1n);

      // Total field in the lth layer: eqs. (1) and (2) in Yang, Appl. Opt., 42 (2003) 1710-1720
      std::complex<double> En = std::pow(c_i, rn)*(2.0*rn + 1.0)/(rn*rn + rn);
      for (int i = 0; i < 3; i++) {
        El[i] = El[i] + En*(cnl_[l][n]*M1o1n[i] - c_i*dnl_[l][n]*N1e1n[i]
                      + c_i*anl_[l][n]*N3e1n[i] - bnl_[l][n]*M3o1n[i]);

        Hl[i] = Hl[i] + En*(-dnl_[l][n]*M1e1n[i] - c_i*cnl_[l][n]*N1o1n[i]
                      + c_i*bnl_[l][n]*N3o1n[i] + anl_[l][n]*M3e1n[i]);
      }
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
    }
   }  // end of MultiLayerMie::fieldInt(...)


  //**********************************************************************************//
  // This function calculates complex electric and magnetic field in the surroundings //
  // and inside (TODO) the particle.                                                  //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send 0 (zero)                    //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
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
  void MultiLayerMie::RunFieldCalculation() {
    // Calculate external scattering coefficients an_ and bn_
    ExtScattCoeffs();
    // Calculate internal scattering coefficients anl_ and bnl_
    IntScattCoeffs();

//    for (int i = 0; i < nmax_; i++) {
//      printf("a[%i] = %g, %g; b[%i] = %g, %g\n", i, an_[i].real(), an_[i].imag(), i, bn_[i].real(), bn_[i].imag());
//    }

    std::vector<double> Pi(nmax_), Tau(nmax_);
    long total_points = coords_[0].size();
    E_.resize(total_points);
    H_.resize(total_points);
    for (auto& f : E_) f.resize(3);
    for (auto& f : H_) f.resize(3);

    for (int point = 0; point < total_points; point++) {
      const double& Xp = coords_[0][point];
      const double& Yp = coords_[1][point];
      const double& Zp = coords_[2][point];

      // Convert to spherical coordinates
      double Rho = std::sqrt(pow2(Xp) + pow2(Yp) + pow2(Zp));

      // If Rho=0 then Theta is undefined. Just set it to zero to avoid problems
      double Theta = (Rho > 0.0) ? std::acos(Zp/Rho) : 0.0;

      // If Xp=Yp=0 then Phi is undefined. Just set it to zero to avoid problems
      double Phi = (Xp != 0.0 || Yp != 0.0) ? std::atan2(Yp, Xp) : 0.0;

      // Avoid convergence problems due to Rho too small
      if (Rho < 1e-10) Rho = 1e-10;

      calcPiTau(std::cos(Theta), Pi, Tau);

      //*******************************************************//
      // external scattering field = incident + scattered      //
      // BH p.92 (4.37), 94 (4.45), 95 (4.50)                  //
      // assume: medium is non-absorbing; refim = 0; Uabs = 0  //
      //*******************************************************//

      // This array contains the fields in spherical coordinates
      std::vector<std::complex<double> > Es(3), Hs(3);

      // Firstly the easiest case: the field outside the particle
      if (Rho > GetSizeParameter()) {
        fieldInt(Rho, Phi, Theta, Pi, Tau, Es, Hs);
//        fieldExt(Rho, Phi, Theta, Pi, Tau, Es, Hs);
        //printf("\nFin E ext: %g,%g,%g   Rho=%g\n", std::abs(Es[0]), std::abs(Es[1]),std::abs(Es[2]), Rho);
      } else {
        fieldInt(Rho, Phi, Theta, Pi, Tau, Es, Hs);
//        printf("\nFin E int: %g,%g,%g   Rho=%g\n", std::abs(Es[0]), std::abs(Es[1]),std::abs(Es[2]), Rho);
      }

      //Now, convert the fields back to cartesian coordinates
      {
        using std::sin;
        using std::cos;
        E_[point][0] = sin(Theta)*cos(Phi)*Es[0] + cos(Theta)*cos(Phi)*Es[1] - sin(Phi)*Es[2];
        E_[point][1] = sin(Theta)*sin(Phi)*Es[0] + cos(Theta)*sin(Phi)*Es[1] + cos(Phi)*Es[2];
        E_[point][2] = cos(Theta)*Es[0] - sin(Theta)*Es[1];

        H_[point][0] = sin(Theta)*cos(Phi)*Hs[0] + cos(Theta)*cos(Phi)*Hs[1] - sin(Phi)*Hs[2];
        H_[point][1] = sin(Theta)*sin(Phi)*Hs[0] + cos(Theta)*sin(Phi)*Hs[1] + cos(Phi)*Hs[2];
        H_[point][2] = cos(Theta)*Hs[0] - sin(Theta)*Hs[1];
      }
      //printf("Cart E: %g,%g,%g   Rho=%g\n", std::abs(Ex), std::abs(Ey),std::abs(Ez), Rho);
    }  // end of for all field coordinates

  }  //  end of MultiLayerMie::RunFieldCalculation()

}  // end of namespace nmie
