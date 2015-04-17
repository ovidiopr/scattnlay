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
#include "bessel.h"
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
  int nMie(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {

    if (x.size() != L || m.size() != L)
        throw std::invalid_argument("Declared number of layers do not fit x and m!");
    if (Theta.size() != nTheta)
        throw std::invalid_argument("Declared number of sample for Theta is not correct!");
    try {
      MultiLayerMie multi_layer_mie;
      multi_layer_mie.SetLayersSize(x);
      multi_layer_mie.SetLayersIndex(m);
      multi_layer_mie.SetAngles(Theta);
      multi_layer_mie.SetPECLayer(pl);
      multi_layer_mie.SetMaxTerms(nmax);

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
  int nMie(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
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
  int nMie(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
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
  int nMie(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
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
  int nField(const unsigned int L, const int pl, const std::vector<double>& x, const std::vector<std::complex<double> >& m, const int nmax, const unsigned int ncoord, const std::vector<double>& Xp_vec, const std::vector<double>& Yp_vec, const std::vector<double>& Zp_vec, std::vector<std::vector<std::complex<double> > >& E, std::vector<std::vector<std::complex<double> > >& H) {
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
      //multi_layer_mie.SetPECLayer(pl); // TODO add PEC layer to field plotting
      multi_layer_mie.SetLayersSize(x);
      multi_layer_mie.SetLayersIndex(m);
      multi_layer_mie.SetFieldCoords({Xp_vec, Yp_vec, Zp_vec});
      multi_layer_mie.RunFieldCalculation();
      E = multi_layer_mie.GetFieldE();
      H = multi_layer_mie.GetFieldH();
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
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    theta_ = angles;
  }


  // ********************************************************************** //
  // Modify size of all layers                                             //
  // ********************************************************************** //
  void MultiLayerMie::SetLayersSize(const std::vector<double>& layer_size) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
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
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    refractive_index_ = index;
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
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    if (layer_position < 0 && layer_position != -1)
      throw std::invalid_argument("Error! Layers are numbered from 0!");
    PEC_layer_position_ = layer_position;
  }


  // ********************************************************************** //
  // Set maximun number of terms to be used                                 //
  // ********************************************************************** //
  void MultiLayerMie::SetMaxTerms(int nmax) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
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
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    size_param_.clear();
    refractive_index_.clear();
  }


  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  //                         Computational core
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //


  // ********************************************************************** //
  // Calculate calcNstop - equation (17)                                    //
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
  void MultiLayerMie::calcNmax(unsigned int first_layer) {
    int ri, riM1;
    const std::vector<double>& x = size_param_;
    const std::vector<std::complex<double> >& m = refractive_index_;
    calcNstop();  // Set initial nmax_ value
    for (unsigned int i = first_layer; i < x.size(); i++) {
      if (static_cast<int>(i) > PEC_layer_position_)  // static_cast used to avoid warning
        ri = round(std::abs(x[i]*m[i]));
      else
        ri = 0;
      nmax_ = std::max(nmax_, ri);
      // first layer is pec, if pec is present
      if ((i > first_layer) && (static_cast<int>(i - 1) > PEC_layer_position_))
        riM1 = round(std::abs(x[i - 1]* m[i]));
      else
        riM1 = 0;
      nmax_ = std::max(nmax_, riM1);
    }
    nmax_ += 15;  // Final nmax_ value
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
  // Calculate an and bn for bulk sphere size x and index m                 //
  // equation (4.56) and (4.57) BH                                          //
  // ********************************************************************** //
  void MultiLayerMie::calc_an_bn_bulk(std::vector<std::complex<double> >& an,
                                      std::vector<std::complex<double> >& bn,
                                      double x, std::complex<double> m) {

    //printf("==========\n m = %g,%g,    x= %g\n", std::real(m), std::imag(m), x);
    std::vector<std::complex<double> > PsiX(nmax_ + 1), ZetaX(nmax_ + 1);
    std::vector<std::complex<double> > PsiMX(nmax_ + 1), ZetaMX(nmax_ + 1);
    // First, calculate the Riccati-Bessel functions
    calcPsiZeta(x, PsiX, ZetaX);
    calcPsiZeta(m*x, PsiMX, ZetaMX);
    std::vector<std::complex<double> > D1X(nmax_ + 1),  D3X(nmax_ + 1);
    std::vector<std::complex<double> > D1MX(nmax_ + 1),  D3MX(nmax_ + 1);
    // Calculate the logarithmic derivatives
    calcD1D3(x, D1X, D3X);
    calcD1D3(m*x, D1MX, D3MX);
    std::vector<std::complex<double> > dPsiX(nmax_ + 1), dZetaX(nmax_ + 1);
    std::vector<std::complex<double> > dPsiMX(nmax_ + 1);
    for (int i = 0; i < nmax_ + 1; ++i) {
      dPsiX[i] = D1X[i]*PsiX[i];
      dPsiMX[i] = D1MX[i]*PsiMX[i];
      //dZetaX[i] = D3X[i]*ZetaX[i];
    }
    bessel::calcZeta(nmax_, x, ZetaX, dZetaX);
    an.resize(nmax_);
    bn.resize(nmax_);
    for (int i = 0; i < nmax_; i++) {
      int n = i+1;
      std::complex<double> Num = m*PsiMX[n]*dPsiX[n] - PsiX[n]*dPsiMX[n];
      std::complex<double> Denom = m*PsiMX[n]*dZetaX[n] - ZetaX[n]*dPsiMX[n];      
      an[i] = Num/Denom;
      Num = PsiMX[n]*dPsiX[n] - m*PsiX[n]*dPsiMX[n];
      Denom = PsiMX[n]*dZetaX[n] - m*ZetaX[n]*dPsiMX[n];      
      bn[i] = Num/Denom;      
    }
    // printf("dPsiX\n");
    // for (auto a: dPsiX) printf("%11.4er%+10.5ei  ",std::real(a), std::imag(a));
    // printf("\ndPsiMX\n");
    // for (auto a: dPsiMX) printf("%11.4er%+10.5ei  ",std::real(a), std::imag(a));
    // printf("\nPsiX\n");
    // for (auto a: PsiX) printf("%11.4er%+10.5ei  ",std::real(a), std::imag(a));
    // printf("\nPsiMX\n");
    // for (auto a: PsiMX) printf("%11.4er%+10.5ei  ",std::real(a), std::imag(a));
    // printf("\nZetaX\n");
    // for (auto a: ZetaX) printf("%11.4er%+10.5ei  ",std::real(a), std::imag(a));
    // printf("\ndZetaX\n");
    // for (auto a: dZetaX) printf("%11.4er%+10.5ei  ",std::real(a), std::imag(a));

    // printf("\nsize param: %g\n", x);


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
      D1[n - 1] = static_cast<double>(n)*zinv - 1.0/(D1[n] + static_cast<double>(n)*zinv);
    }

    if (std::abs(D1[0]) > 100000.0)
      throw std::invalid_argument("Unstable D1! Please, try to change input parameters!\n");

    // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
    PsiZeta_[0] = 0.5*(1.0 - std::complex<double>(std::cos(2.0*z.real()), std::sin(2.0*z.real()))
                      *std::exp(-2.0*z.imag()));
    D3[0] = std::complex<double>(0.0, 1.0);
    for (int n = 1; n <= nmax_; n++) {
      PsiZeta_[n] = PsiZeta_[n - 1]*(static_cast<double>(n)*zinv - D1[n - 1])
                                   *(static_cast<double>(n)*zinv - D3[n - 1]);
      D3[n] = D1[n] + std::complex<double>(0.0, 1.0)/PsiZeta_[n];
    }
  }


  //**********************************************************************************//
  // This function calculates the Riccati-Bessel functions (Psi and Zeta) for a       //
  // complex argument (z).                                                            //
  // Equations (20a) - (21b)                                                          //
  //                                                                                  //
  // Input parameters:                                                                //
  //   z: Complex argument to evaluate Psi and Zeta                                   //
  //   nmax: Maximum number of terms to calculate Psi and Zeta                        //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Psi, Zeta: Riccati-Bessel functions                                            //
  //**********************************************************************************//
  void MultiLayerMie::calcPsiZeta(std::complex<double> z,
                                  std::vector<std::complex<double> >& Psi,
                                  std::vector<std::complex<double> >& Zeta) {

    std::complex<double> c_i(0.0, 1.0);
    std::vector<std::complex<double> > D1(nmax_ + 1), D3(nmax_ + 1);

    // First, calculate the logarithmic derivatives
    calcD1D3(z, D1, D3);

    // Now, use the upward recurrence to calculate Psi and Zeta - equations (20a) - (21b)
    Psi[0] = std::sin(z);
    Zeta[0] = std::sin(z) - c_i*std::cos(z);
    for (int n = 1; n <= nmax_; n++) {
      Psi[n]  =  Psi[n - 1]*(static_cast<double>(n)/z - D1[n - 1]);
      Zeta[n] = Zeta[n - 1]*(static_cast<double>(n)/z - D3[n - 1]);
    }
  }


  //**********************************************************************************//
  // This function calculates the spherical Bessel (jn) and Hankel (h1n) functions    //
  // and their derivatives for a given complex value z. See pag. 87 B&H.              //
  //                                                                                  //
  // Input parameters:                                                                //
  //   z: Complex argument to evaluate jn and h1n                                     //
  //   nmax_: Maximum number of terms to calculate jn and h1n                         //
  //                                                                                  //
  // Output parameters:                                                               //
  //   jn, h1n: Spherical Bessel and Hankel functions                                 //
  //   jnp, h1np: Derivatives of the spherical Bessel and Hankel functions            //
  //                                                                                  //
  // What we actually calculate are the Ricatti-Bessel fucntions and then simply      //
  // evaluate the spherical Bessel and Hankel functions and their derivatives         //
  // using the relations:                                                             //
  //                                                                                  //
  //     j[n]   = Psi[n]/z                                                            //
  //     j'[n]  = j[n-1]-(n+1)*jn[n])/z                                               //
  //     h1[n]  = Zeta[n]/z                                                           //
  //     h1'[n] = h1[n-1]-(n+1)*h1n[n]/z                                              //
  //                                                                                  //
  //**********************************************************************************//
  void MultiLayerMie::sbesjh(std::complex<double> z,
                             std::vector<std::complex<double> >& jn, std::vector<std::complex<double> >& jnp,
                             std::vector<std::complex<double> >& h1n, std::vector<std::complex<double> >& h1np) {

    std::vector<std::complex<double> > Psi(nmax_ + 1), Zeta(nmax_ + 1);

    // First, calculate the Riccati-Bessel functions
    calcPsiZeta(z, Psi, Zeta);

    // Now, calculate Spherical Bessel and Hankel functions and their derivatives
    for (int n = 0; n <= nmax_; n++) {
      jn[n] = Psi[n]/z;
      h1n[n] = Zeta[n]/z;

      if (n == 0) {
        jnp[0] = -Psi[1]/z - jn[0]/z;
        h1np[0] = -Zeta[1]/z - h1n[0]/z;
      } else {
        jnp[n] = jn[n - 1] - static_cast<double>(n + 1)*jn[n]/z;
        h1np[n] = h1n[n - 1] - static_cast<double>(n + 1)*h1n[n]/z;
      }
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

    int i;
    //****************************************************//
    // Equations (26a) - (26c)                            //
    //****************************************************//
    // Initialize Pi and Tau
    Pi[0] = 1.0;  // n=1
    Tau[0] = costheta;
    // Calculate the actual values
    if (nmax_ > 1) {
      Pi[1] = 3*costheta*Pi[0]; //n=2
      Tau[1] = 2*costheta*Pi[1] - 3*Pi[0];
      for (i = 2; i < nmax_; i++) { //n=[3..nmax_]
        Pi[i] = ((i + i + 1)*costheta*Pi[i - 1] - (i + 1)*Pi[i - 2])/i;
        Tau[i] = (i + 1)*costheta*Pi[i] - (i + 2)*Pi[i - 1];
      }
    }
  }  // end of MultiLayerMie::calcPiTau(...)


  //**********************************************************************************//
  // This function calculates vector spherical harmonics (eq. 4.50, p. 95 BH),        //
  // required to calculate the near-field parameters.                                 //
  //                                                                                  //
  // Input parameters:                                                                //
  //   Rho: Radial distance                                                           //
  //   Phi: Azimuthal angle                                                           //
  //   Theta: Polar angle                                                             //
  //   rn: Either the spherical Ricatti-Bessel function of first or third kind        //
  //   Dn: Logarithmic derivative of rn                                               //
  //   Pi, Tau: Angular functions Pi and Tau                                          //
  //   n: Order of vector spherical harmonics                                         //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Mo1n, Me1n, No1n, Ne1n: Complex vector spherical harmonics                     //
  //**********************************************************************************//
  void MultiLayerMie::calcSpherHarm(const double Rho, const double Theta, const double Phi,
                                    const std::complex<double>& rn, const std::complex<double>& Dn,
                                    const double& Pi, const double& Tau, const double& n,
                                    std::vector<std::complex<double> >& Mo1n, std::vector<std::complex<double> >& Me1n, 
                                    std::vector<std::complex<double> >& No1n, std::vector<std::complex<double> >& Ne1n) {

    // using eq 4.50 in BH
    std::complex<double> c_zero(0.0, 0.0);

    using std::sin;
    using std::cos;
    Mo1n[0] = c_zero;
    Mo1n[1] = cos(Phi)*Pi*rn/Rho;
    Mo1n[2] = -sin(Phi)*Tau*rn/Rho;
    Me1n[0] = c_zero;
    Me1n[1] = -sin(Phi)*Pi*rn/Rho;
    Me1n[2] = -cos(Phi)*Tau*rn/Rho;
    No1n[0] = sin(Phi)*(n*n + n)*sin(Theta)*Pi*rn/Rho/Rho;
    No1n[1] = sin(Phi)*Tau*Dn*rn/Rho;
    No1n[2] = cos(Phi)*Pi*Dn*rn/Rho;
    Ne1n[0] = cos(Phi)*(n*n + n)*sin(Theta)*Pi*rn/Rho/Rho;
    Ne1n[1] = cos(Phi)*Tau*Dn*rn/Rho;
    Ne1n[2] = -sin(Phi)*Pi*Dn*rn/Rho;
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
  void MultiLayerMie::ScattCoeffs() {

    isScaCoeffsCalc_ = false;

    const std::vector<double>& x = size_param_;
    const std::vector<std::complex<double> >& m = refractive_index_;
    const int& pl = PEC_layer_position_;
    const int L = refractive_index_.size();

    //************************************************************************//
    // Calculate the index of the first layer. It can be either 0 (default)   //
    // or the index of the outermost PEC layer. In the latter case all layers //
    // below the PEC are discarded.                                           //
    // ***********************************************************************//
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

    std::vector<std::complex<double> > PsiXL(nmax_ + 1), ZetaXL(nmax_ + 1);

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
      //Calculate D1 and D3 for z1 and z2 in the layers fl + 1..L   //
      //************************************************************//
      z1 = x[l]*m[l];
      z2 = x[l - 1]*m[l];
      //Calculate D1 and D3 for z1
      calcD1D3(z1, D1_mlxl, D3_mlxl);
      //Calculate D1 and D3 for z2
      calcD1D3(z2, D1_mlxlM1, D3_mlxlM1);

      //*************************************************//
      //Calculate Q, Ha and Hb in the layers fl + 1..L   //
      //*************************************************//
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
    //Calculate Psi and Zeta for XL         //
    //**************************************//
    // Calculate PsiXL and ZetaXL
    calcPsiZeta(x[L - 1], PsiXL, ZetaXL);

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
    isScaCoeffsCalc_ = true;
  }  // end of MultiLayerMie::ScattCoeffs(...)


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
    if (size_param_.size() != refractive_index_.size())
      throw std::invalid_argument("Each size parameter should have only one index!");
    if (size_param_.size() == 0)
      throw std::invalid_argument("Initialize model first!");

    const std::vector<double>& x = size_param_;

    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;

    // Calculate scattering coefficients
    ScattCoeffs();

    if (!isScaCoeffsCalc_) // TODO seems to be unreachable
      throw std::invalid_argument("Calculation of scattering coefficients failed!");

    // Initialize the scattering parameters
    Qext_ = 0.0;
    Qsca_ = 0.0;
    Qabs_ = 0.0;
    Qbk_ = 0.0;
    Qpr_ = 0.0;
    asymmetry_factor_ = 0.0;
    albedo_ = 0.0;

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
      Qext_ += (n + n + 1.0)*(an_[i].real() + bn_[i].real());
      // Equation (28)
      Qsca_ += (n + n + 1.0)*(an_[i].real()*an_[i].real() + an_[i].imag()*an_[i].imag()
                  + bn_[i].real()*bn_[i].real() + bn_[i].imag()*bn_[i].imag());
      // Equation (29)
      Qpr_ += ((n*(n + 2)/(n + 1))*((an_[i]*std::conj(an_[n]) + bn_[i]*std::conj(bn_[n])).real())
               + ((double)(n + n + 1)/(n*(n + 1)))*(an_[i]*std::conj(bn_[i])).real());
      // Equation (33)
      Qbktmp += (double)(n + n + 1)*(1 - 2*(n % 2))*(an_[i]- bn_[i]);
      // Calculate the scattering amplitudes (S1 and S2)    //
      // Equations (25a) - (25b)                            //
      for (unsigned int t = 0; t < theta_.size(); t++) {
        calcPiTau(std::cos(theta_[t]), Pi, Tau);

        S1_[t] += calc_S1(n, an_[i], bn_[i], Pi[i], Tau[i]);
        S2_[t] += calc_S2(n, an_[i], bn_[i], Pi[i], Tau[i]);
      }
    }
    double x2 = pow2(x.back());
    Qext_ = 2.0*(Qext_)/x2;                                 // Equation (27)
    Qsca_ = 2.0*(Qsca_)/x2;                                 // Equation (28)
    Qpr_ = Qext_ - 4.0*(Qpr_)/x2;                           // Equation (29)
    Qabs_ = Qext_ - Qsca_;                                  // Equation (30)
    albedo_ = Qsca_/Qext_;                                  // Equation (31)
    asymmetry_factor_ = (Qext_ - Qpr_)/Qsca_;               // Equation (32)
    Qbk_ = (Qbktmp.real()*Qbktmp.real() + Qbktmp.imag()*Qbktmp.imag())/x2;    // Equation (33)

    isMieCalculated_ = true;
  }


  //**********************************************************************************//
  // This function calculates the expansion coefficients inside the particle,         //
  // required to calculate the near-field parameters.                                 //
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
  //   aln, bln, cln, dln: Complex scattering amplitudes inside the particle          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  void MultiLayerMie::ExpanCoeffs() {
    if (!isScaCoeffsCalc_)
      throw std::invalid_argument("(ExpanCoeffs) You should calculate external coefficients first!");

    isExpCoeffsCalc_ = false;

    std::complex<double> c_one(1.0, 0.0), c_zero(0.0, 0.0);

    const int L = refractive_index_.size();

    aln_.resize(L + 1);
    bln_.resize(L + 1);
    cln_.resize(L + 1);
    dln_.resize(L + 1);
    for (int l = 0; l <= L; l++) {
      aln_[l].resize(nmax_);
      bln_[l].resize(nmax_);
      cln_[l].resize(nmax_);
      dln_[l].resize(nmax_);
    }

    // Yang, paragraph under eq. A3
    // a^(L + 1)_n = a_n, d^(L + 1) = 1 ...
    for (int n = 0; n < nmax_; n++) {
      aln_[L][n] = an_[n];
      bln_[L][n] = bn_[n];
      cln_[L][n] = c_one;
      dln_[L][n] = c_one;

      printf("aln_[%02i, %02i] = %g,%g; bln_[%02i, %02i] = %g,%g; cln_[%02i, %02i] = %g,%g; dln_[%02i, %02i] = %g,%g\n", L, n, std::real(aln_[L][n]), std::imag(aln_[L][n]), L, n, std::real(bln_[L][n]), std::imag(bln_[L][n]), L, n, std::real(cln_[L][n]), std::imag(cln_[L][n]), L, n, std::real(dln_[L][n]), std::imag(dln_[L][n]));
    }

    std::vector<std::complex<double> > D1z(nmax_ + 1), D1z1(nmax_ + 1), D3z(nmax_ + 1), D3z1(nmax_ + 1);
    std::vector<std::complex<double> > Psiz(nmax_ + 1), Psiz1(nmax_ + 1), Zetaz(nmax_ + 1), Zetaz1(nmax_ + 1);
    std::complex<double> denomZeta, denomPsi, T1, T2, T3, T4;

    auto& m = refractive_index_;
    std::vector< std::complex<double> > m1(L);

    for (int l = 0; l < L - 1; l++) m1[l] = m[l + 1];
    m1[L - 1] = std::complex<double> (1.0, 0.0);

    std::complex<double> z, z1;
    for (int l = L - 1; l >= 0; l--) {
      z = size_param_[l]*m[l];
      z1 = size_param_[l]*m1[l];

      calcD1D3(z, D1z, D3z);
      calcD1D3(z1, D1z1, D3z1);
      calcPsiZeta(z, Psiz, Zetaz);
      calcPsiZeta(z1, Psiz1, Zetaz1);

      for (int n = 0; n < nmax_; n++) {
        int n1 = n + 1;

        denomZeta = Zetaz[n1]*(D1z[n1] - D3z[n1]);
        denomPsi  =  Psiz[n1]*(D1z[n1] - D3z[n1]);

        T1 = aln_[l + 1][n]*Zetaz1[n1] - dln_[l + 1][n]*Psiz1[n1];
        T2 = (bln_[l + 1][n]*Zetaz1[n1] - cln_[l + 1][n]*Psiz1[n1])*m[l]/m1[l];

        T3 = (D1z1[n1]*dln_[l + 1][n]*Psiz1[n1] - D3z1[n1]*aln_[l + 1][n]*Zetaz1[n1])*m[l]/m1[l];
        T4 = D1z1[n1]*cln_[l + 1][n]*Psiz1[n1] - D3z1[n1]*bln_[l + 1][n]*Zetaz1[n1];

        // aln
        aln_[l][n] = (D1z[n1]*T1 + T3)/denomZeta;
        // bln
        bln_[l][n] = (D1z[n1]*T2 + T4)/denomZeta;
        // cln
        cln_[l][n] = (D3z[n1]*T2 + T4)/denomPsi;
        // dln
        dln_[l][n] = (D3z[n1]*T1 + T3)/denomPsi;

        printf("aln_[%02i, %02i] = %g,%g; bln_[%02i, %02i] = %g,%g; cln_[%02i, %02i] = %g,%g; dln_[%02i, %02i] = %g,%g\n", l, n, real(aln_[l][n]), imag(aln_[l][n]), l, n, real(bln_[l][n]), imag(bln_[l][n]), l, n, real(cln_[l][n]), imag(cln_[l][n]), l, n, real(dln_[l][n]), imag(dln_[l][n]));
      }  // end of all n
    }  // end of all l

    // Check the result and change  aln_[0][n] and aln_[0][n] for exact zero
    for (int n = 0; n < nmax_; ++n) {
//      printf("n=%d, aln_=%g,%g,   bln_=%g,%g \n", n, real(aln_[0][n]), imag(aln_[0][n]),
//         real(bln_[0][n]), imag(bln_[0][n]));
      if (std::abs(aln_[0][n]) < 1e-10) aln_[0][n] = 0.0;
      else throw std::invalid_argument("Unstable calculation of aln_[0][n]!");
      if (std::abs(bln_[0][n]) < 1e-10) bln_[0][n] = 0.0;
      else throw std::invalid_argument("Unstable calculation of bln_[0][n]!");
    }

    isExpCoeffsCalc_ = true;
  }  // end of   void MultiLayerMie::ExpanCoeffs()


  //**********************************************************************************//
  // This function calculates the expansion coefficients inside the particle,         //
  // required to calculate the near-field parameters.                                 //
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
  //   aln, bln, cln, dln: Complex scattering amplitudes inside the particle          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  void MultiLayerMie::ExpanCoeffsV2() {
    if (!isScaCoeffsCalc_)
      throw std::invalid_argument("(ExpanCoeffs) You should calculate external coefficients first!");

    isExpCoeffsCalc_ = false;

    std::complex<double> c_one(1.0, 0.0), c_zero(0.0, 0.0);

    const int L = refractive_index_.size();

    aln_.resize(L + 1);
    bln_.resize(L + 1);
    cln_.resize(L + 1);
    dln_.resize(L + 1);
    for (int l = 0; l <= L; l++) {
      aln_[l].resize(nmax_);
      bln_[l].resize(nmax_);
      cln_[l].resize(nmax_);
      dln_[l].resize(nmax_);
    }

    // Yang, paragraph under eq. A3
    // a^(L + 1)_n = a_n, d^(L + 1) = 1 ...
    for (int n = 0; n < nmax_; n++) {
      aln_[L][n] = an_[n];
      bln_[L][n] = bn_[n];
      cln_[L][n] = c_one;
      dln_[L][n] = c_one;

      printf("aln_[%02i, %02i] = %g,%g; bln_[%02i, %02i] = %g,%g; cln_[%02i, %02i] = %g,%g; dln_[%02i, %02i] = %g,%g\n", L, n, std::real(aln_[L][n]), std::imag(aln_[L][n]), L, n, std::real(bln_[L][n]), std::imag(bln_[L][n]), L, n, std::real(cln_[L][n]), std::imag(cln_[L][n]), L, n, real(dln_[L][n]), std::imag(dln_[L][n]));
    }

    std::vector<std::complex<double> > D1z(nmax_ + 1), D1z1(nmax_ + 1), D3z(nmax_ + 1), D3z1(nmax_ + 1);
    std::vector<std::complex<double> > Psiz(nmax_ + 1), Psiz1(nmax_ + 1), Zetaz(nmax_ + 1), Zetaz1(nmax_ + 1);
    std::complex<double> denomZeta, denomPsi, T1, T2, T3, T4;

    std::vector<std::vector<std::complex<double> > > a(2);
    a[0].resize(3);
    a[1].resize(3);

    auto& m = refractive_index_;
    std::vector< std::complex<double> > m1(L);

    for (int l = 0; l < L - 1; l++) m1[l] = m[l + 1];
    m1[L - 1] = std::complex<double> (1.0, 0.0);

    std::complex<double> z, z1;
    for (int l = L - 1; l >= 0; l--) {
      z = size_param_[l]*m[l];
      z1 = size_param_[l]*m1[l];

      calcD1D3(z, D1z, D3z);
      calcD1D3(z1, D1z1, D3z1);
      calcPsiZeta(z, Psiz, Zetaz);
      calcPsiZeta(z1, Psiz1, Zetaz1);

      for (int n = 0; n < nmax_; n++) {
        int n1 = n + 1;

        a[0][0] = m1[l]*D3z[n1]*Zetaz[n1];
        a[0][1] = -m1[l]*D1z[n1]*Psiz[n1];
        a[0][2] = aln_[l + 1][n]*m[l]*D3z1[n1]*Zetaz1[n1];
        a[0][2] -= dln_[l + 1][n]*m[l]*D1z1[n1]*Psiz1[n1];

        a[1][0] = Zetaz[n1];
        a[1][1] = -Psiz[n1];
        a[1][2] = aln_[l + 1][n]*Zetaz1[n1] - dln_[l + 1][n]*Psiz1[n1];

        // aln
        aln_[l][n] = (a[0][2]*a[1][1] - a[0][1]*a[1][2])/(a[0][0]*a[1][1] - a[0][1]*a[1][0]);
        // dln
        dln_[l][n] = (a[0][2]*a[1][0] - a[0][0]*a[1][2])/(a[0][1]*a[1][0] - a[0][0]*a[0][1]);

        /*for (int i = 0; i < 2; i++) {
          for (int j = 0; j < 3; j++) {
            printf("a[%i, %i] = %g,%g ", i, j, real(a[i][j]), imag(a[i][j]));
          }
          printf("\n");
        }
        printf("aln_[%i, %i] = %g,%g; dln_[%i, %i] = %g,%g\n\n", l, n, real(aln_[l][n]), imag(aln_[l][n]), l, n, real(dln_[l][n]), imag(dln_[l][n]));*/

        a[0][0] = D3z[n1]*Zetaz[n1];
        a[0][1] = -D1z[n1]*Psiz[n1];
        a[0][2] = bln_[l + 1][n]*D3z1[n1]*Zetaz1[n1];
        a[0][2] -= cln_[l + 1][n]*D1z1[n1]*Psiz1[n1];

        a[1][0] = m1[l]*Zetaz[n1];
        a[1][1] = -m1[l]*Psiz[n1];
        a[1][2] = bln_[l + 1][n]*m[l]*Zetaz1[n1] - cln_[l + 1][n]*m[l]*Psiz1[n1];

        // bln
        bln_[l][n] = (a[0][2]*a[1][1] - a[0][1]*a[1][2])/(a[0][0]*a[1][1] - a[0][1]*a[1][0]);
        // cln
        cln_[l][n] = (a[0][2]*a[1][0] - a[0][0]*a[1][2])/(a[0][1]*a[1][0] - a[0][0]*a[0][1]);

        printf("aln_[%02i, %02i] = %g,%g; bln_[%02i, %02i] = %g,%g; cln_[%02i, %02i] = %g,%g; dln_[%02i, %02i] = %g,%g\n", l, n, real(aln_[l][n]), imag(aln_[l][n]), l, n, real(bln_[l][n]), imag(bln_[l][n]), l, n, real(cln_[l][n]), imag(cln_[l][n]), l, n, real(dln_[l][n]), imag(dln_[l][n]));
      }  // end of all n
    }  // end of all l

    // Check the result and change  aln_[0][n] and aln_[0][n] for exact zero
    for (int n = 0; n < nmax_; ++n) {
//      printf("n=%d, aln_=%g,%g,   bln_=%g,%g \n", n, real(aln_[0][n]), imag(aln_[0][n]),
//         real(bln_[0][n]), imag(bln_[0][n]));
      if (std::abs(aln_[0][n]) < 1e-1) aln_[0][n] = 0.0;
      else throw std::invalid_argument("Unstable calculation of aln_[0][n]!");
      if (std::abs(bln_[0][n]) < 1e-1) bln_[0][n] = 0.0;
      else throw std::invalid_argument("Unstable calculation of bln_[0][n]!");
    }

    isExpCoeffsCalc_ = true;
  }  // end of   void MultiLayerMie::ExpanCoeffs()


  // ********************************************************************** //
  // external scattering field = incident + scattered                       //
  // BH p.92 (4.37), 94 (4.45), 95 (4.50)                                   //
  // assume: medium is non-absorbing; refim = 0; Uabs = 0                   //
  // ********************************************************************** //
  void MultiLayerMie::fieldExt(const double Rho, const double Theta, const double Phi,
                               std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H)  {

    std::complex<double> c_zero(0.0, 0.0), c_i(0.0, 1.0), c_one(1.0, 0.0);
    std::vector<std::complex<double> > ipow = {c_one, c_i, -c_one, -c_i}; // Vector containing precomputed integer powers of i to avoid computation
    std::vector<std::complex<double> > M3o1n(3), M3e1n(3), N3o1n(3), N3e1n(3);
    std::vector<std::complex<double> > Ei(3, c_zero), Hi(3, c_zero), Es(3, c_zero), Hs(3, c_zero);
    std::vector<std::complex<double> > Psi(nmax_ + 1), D1n(nmax_ + 1), Zeta(nmax_ + 1), D3n(nmax_ + 1);
    std::vector<double> Pi(nmax_), Tau(nmax_);

    // Avoid calculation inside the particle
    if (Rho < size_param_.back()) {
       for (int i = 0; i < 3; i++) {
        E[i] = c_zero;
        H[i] = c_zero;
      }
      return;
    }

    calcD1D3(Rho, D1n, D3n);
    calcPsiZeta(Rho, Psi, Zeta);

    // Calculate spherical Bessel and Hankel functions
    //sbesjh(Rho, jn, jnp, h1n, h1np);

    // Calculate angular functions Pi and Tau
    calcPiTau(std::cos(Theta), Pi, Tau);

    for (int n = 0; n < nmax_; n++) {
      int n1 = n + 1;
      double rn = static_cast<double>(n1);

      // using BH 4.12 and 4.50
      calcSpherHarm(Rho, Theta, Phi, Zeta[n1], D3n[n1], Pi[n], Tau[n], rn, M3o1n, M3e1n, N3o1n, N3e1n);

      // scattered field: BH p.94 (4.45)
      std::complex<double> En = ipow[n1 % 4]*(rn + rn + 1.0)/(rn*rn + rn);
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
    }
   }  // end of MultiLayerMie::fieldExt(...)


  //**********************************************************************************//
  // This function calculates the electric (E) and magnetic (H) fields inside and     //
  // around the particle.                                                             //
  //                                                                                  //
  // Input parameters (coordinates of the point):                                     //
  //   Rho: Radial distance                                                           //
  //   Phi: Azimuthal angle                                                           //
  //   Theta: Polar angle                                                             //
  //                                                                                  //
  // Output parameters:                                                               //
  //   E, H: Complex electric and magnetic fields                                     //
  //**********************************************************************************//
  void MultiLayerMie::calcField(const double Rho, const double Theta, const double Phi,
                                std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H)  {

    std::complex<double> c_zero(0.0, 0.0), c_i(0.0, 1.0), c_one(1.0, 0.0);
    std::vector<std::complex<double> > ipow = {c_one, c_i, -c_one, -c_i}; // Vector containing precomputed integer powers of i to avoid computation
    std::vector<std::complex<double> > M3o1n(3), M3e1n(3), N3o1n(3), N3e1n(3);
    std::vector<std::complex<double> > M1o1n(3), M1e1n(3), N1o1n(3), N1e1n(3);
    std::vector<std::complex<double> > Psi(nmax_ + 1), D1n(nmax_ + 1), Zeta(nmax_ + 1), D3n(nmax_ + 1);
    std::vector<double> Pi(nmax_), Tau(nmax_);

    std::vector<std::complex<double> > Ei(3), Hi(3);

    int l = 0;  // Layer number
    std::complex<double> ml;

    // Initialize E and H
    for (int i = 0; i < 3; i++) {
      E[i] = c_zero;
      H[i] = c_zero;
    }

    if (Rho > size_param_.back()) {
      l = size_param_.size();
      ml = c_one;
    } else {
      for (int i = size_param_.size() - 1; i >= 0 ; i--) {
        if (Rho <= size_param_[i]) {
          l = i;
      break;
        }
      }
      ml = refractive_index_[l];
    }


    calcD1D3(Rho, D1n, D3n);
    calcPsiZeta(Rho, Psi, Zeta);

    // Calculate spherical Bessel and Hankel functions and their derivatives
    //sbesjh(Rho*ml, jn, jnp, h1n, h1np);
    //sbesjh(2.0*PI_*Rho*ml, jn, jnp, h1n, h1np);
    //printf("2.0*PI*Rho*ml = %10.5er%+10.5ei\n",std::real(2.0*PI_*Rho*ml), std::imag(2.0*PI_*Rho*ml));

    // Calculate angular functions Pi and Tau
    calcPiTau(std::cos(Theta), Pi, Tau);
    //printf("Thetd = %g, cos(Theta) = %g\n", Theta, std::cos(Theta));
    //printf("jn:\n");
    //for (auto p : jn) printf("%+11.4er%+11.4ei\n",p.real(), p.imag());
    //printf("Pi:\n");
    //for (auto p : Pi) printf("%11.4e\n",p);
    //printf("Tau:\n");
    //for (auto p : Tau) printf("%11.4e\n",p);

    for (int n = nmax_ - 2; n >= 0; n--) {
      int n1 = n + 1;
      double rn = static_cast<double>(n1);

      // using BH 4.12 and 4.50
      calcSpherHarm(Rho, Theta, Phi, Psi[n1], D1n[n1], Pi[n], Tau[n], rn, M1o1n, M1e1n, N1o1n, N1e1n);
      calcSpherHarm(Rho, Theta, Phi, Zeta[n1], D3n[n1], Pi[n], Tau[n], rn, M3o1n, M3e1n, N3o1n, N3e1n);
//      auto deriv1 = -rn*jn[n1]+Rho*jn[n1-1];
//      auto deriv2 = Rho*jnp[n1] + jn[n1];
//      printf("n=%d   deriv1: %+11.4e   deriv2: %+11.4ei\n",n1, deriv1.real(), deriv2.real());
      // printf("N1e1n[%d]: ", n1);
      // for (auto p : N1e1n) printf("%+11.4er%+11.4ei\t",p.real(), p.imag());
      // printf("\n");

      // Total field in the lth layer: eqs. (1) and (2) in Yang, Appl. Opt., 42 (2003) 1710-1720
      std::complex<double> En = ipow[n1 % 4]*(rn + rn + 1.0)/(rn*rn + rn);
      for (int i = 0; i < 3; i++) {
        // electric field E [V m - 1] = EF*E0
        E[i] += En*(cln_[l][n]*M1o1n[i] - c_i*dln_[l][n]*N1e1n[i]
              + c_i*aln_[l][n]*N3e1n[i] -     bln_[l][n]*M3o1n[i]);

        H[i] += En*(-dln_[l][n]*M1e1n[i] - c_i*cln_[l][n]*N1o1n[i]
              +  c_i*bln_[l][n]*N3o1n[i] +     aln_[l][n]*M3e1n[i]);
        Ei[i] += En*(M1o1n[i] - c_i*N1e1n[i]);
        Hi[i] += En*(-M1e1n[i] - c_i*N1o1n[i]);

      }
    }  // end of for all n

    //printf("rho = %11.4e; phi = %11.4eº; theta = %11.4eº; x[%i] = %11.4e; m[%i] = %11.4er%+10.5ei\n", Rho, Phi*180./PI_, Theta*180./PI_, l, size_param_[l], l, std::real(ml), std::imag(ml));
    // magnetic field
    double hffact = 1.0/(cc_*mu_);
    for (int i = 0; i < 3; i++) {
      H[i] = hffact*H[i];
      Hi[i] *= hffact;
      //printf("E[%i] = %10.5er%+10.5ei; Ei[%i] = %10.5er%+10.5ei; H[%i] = %10.5er%+10.5ei; Hi[%i] = %10.5er%+10.5ei\n", i, std::real(E[i]), std::imag(E[i]), i, std::real(Ei[i]), std::imag(Ei[i]), i, std::real(H[i]), std::imag(H[i]), i, std::real(Hi[i]), std::imag(Hi[i]));
    }
   }  // end of MultiLayerMie::calcField(...)


  //**********************************************************************************//
  // This function calculates complex electric and magnetic field in the surroundings //
  // and inside the particle.                                                         //
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
    double Rho, Theta, Phi;

    // Calculate scattering coefficients an_ and bn_
    ScattCoeffs();

    // std::vector<std::complex<double> > an1(nmax_), bn1(nmax_);
    // calc_an_bn_bulk(an1, bn1, size_param_.back(), refractive_index_.back());
    // for (int n = 0; n < nmax_; n++) {
    //   printf("an_[%i] = %11.4er%+10.5ei;  an_bulk_[%i] = %11.4er%+10.5ei\n", n, std::real(an_[n]), std::imag(an_[n]), n, std::real(an1[n]), std::imag(an1[n]));
    //   printf("bn_[%i] = %11.4er%+10.5ei;  bn_bulk_[%i] = %11.4er%+10.5ei\n", n, std::real(bn_[n]), std::imag(bn_[n]), n, std::real(bn1[n]), std::imag(bn1[n]));
    // }

    // Calculate expansion coefficients aln_,  bln_, cln_, and dln_
    ExpanCoeffs();
    // for (int i = 0; i < nmax_; ++i) {
    //   printf("cln_[%i] = %11.4er%+10.5ei;  dln_[%i] = %11.4er%+10.5ei\n", i, std::real(cln_[0][i]), std::imag(cln_[0][i]),
    //          i, std::real(dln_[0][i]), std::imag(dln_[0][i]));
    // }


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
      Rho = std::sqrt(pow2(Xp) + pow2(Yp) + pow2(Zp));

      // If Rho=0 then Theta is undefined. Just set it to zero to avoid problems
      Theta = (Rho > 0.0) ? std::acos(Zp/Rho) : 0.0;

      // If Xp=Yp=0 then Phi is undefined. Just set it to zero to avoid problems
      if (Xp == 0.0)
        Phi = (Yp != 0.0) ? std::asin(Yp/std::sqrt(pow2(Xp) + pow2(Yp))) : 0.0;
      else
        Phi = std::acos(Xp/std::sqrt(pow2(Xp) + pow2(Yp)));

      // Avoid convergence problems due to Rho too small
      if (Rho < 1e-5) Rho = 1e-5;

      //printf("X = %g; Y = %g; Z = %g; pho = %g; phi = %g; theta = %g\n", Xp, Yp, Zp, Rho, Phi*180./PI_, Theta*180./PI_);

      //*******************************************************//
      // external scattering field = incident + scattered      //
      // BH p.92 (4.37), 94 (4.45), 95 (4.50)                  //
      // assume: medium is non-absorbing; refim = 0; Uabs = 0  //
      //*******************************************************//

      // This array contains the fields in spherical coordinates
      std::vector<std::complex<double> > Es(3), Hs(3);

      // Firstly the easiest case: the field outside the particle
      //      if (Rho >= GetSizeParameter()) {
      //        fieldExt(Rho, Theta, Phi, Es, Hs);
      // } else {
      calcField(Rho, Theta, Phi, Es, Hs);  //Should work fine both: inside and outside the particle
      //}

      { //Now, convert the fields back to cartesian coordinates
        using std::sin;
        using std::cos;
        E_[point][0] = sin(Theta)*cos(Phi)*Es[0] + cos(Theta)*cos(Phi)*Es[1] - sin(Phi)*Es[2];
        E_[point][1] = sin(Theta)*sin(Phi)*Es[0] + cos(Theta)*sin(Phi)*Es[1] + cos(Phi)*Es[2];
        E_[point][2] = cos(Theta)*Es[0] - sin(Theta)*Es[1];

        H_[point][0] = sin(Theta)*cos(Phi)*Hs[0] + cos(Theta)*cos(Phi)*Hs[1] - sin(Phi)*Hs[2];
        H_[point][1] = sin(Theta)*sin(Phi)*Hs[0] + cos(Theta)*sin(Phi)*Hs[1] + cos(Phi)*Hs[2];
        H_[point][2] = cos(Theta)*Hs[0] - sin(Theta)*Hs[1];
      }
    }  // end of for all field coordinates
  }  //  end of MultiLayerMie::RunFieldCalculation()
}  // end of namespace nmie
