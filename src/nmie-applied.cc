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
#include "nmie-applied.hpp"
#include <array>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace nmie {  
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
  int nMieApplied(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {

    if (x.size() != L || m.size() != L)
        throw std::invalid_argument("Declared number of layers do not fit x and m!");
    if (Theta.size() != nTheta)
        throw std::invalid_argument("Declared number of sample for Theta is not correct!");
    try {
      MultiLayerMieApplied<> ml_mie;
      ml_mie.SetLayersSize(x);
      ml_mie.SetLayersIndex(m);
      ml_mie.SetAngles(Theta);
      ml_mie.SetPECLayer(pl);
      ml_mie.SetMaxTerms(nmax);

      ml_mie.RunMieCalculation();

      *Qext = ml_mie.GetQext();
      *Qsca = ml_mie.GetQsca();
      *Qabs = ml_mie.GetQabs();
      *Qbk = ml_mie.GetQbk();
      *Qpr = ml_mie.GetQpr();
      *g = ml_mie.GetAsymmetryFactor();
      *Albedo = ml_mie.GetAlbedo();
      S1 = ml_mie.GetS1();
      S2 = ml_mie.GetS2();

      return ml_mie.GetMaxTerms();
    } catch(const std::invalid_argument& ia) {
      // Will catch if  ml_mie fails or other errors.
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      throw std::invalid_argument(ia);
      return -1;
    }
    return 0;
  }
  int nMieApplied(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
    return nmie::nMieApplied(L, -1, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }
  int nMieApplied(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
    return nmie::nMieApplied(L, pl, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }
  int nMieApplied(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
    return nmie::nMieApplied(L, -1, x, m, nTheta, Theta, nmax, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::GetFailed() {
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
  void MultiLayerMieApplied::AddTargetLayer(double width, std::complex<double> layer_index) {
    MarkUncalculated();
    if (width <= 0)
      throw std::invalid_argument("Layer width should be positive!");
    target_width_.push_back(width);
    target_index_.push_back(layer_index);
  }  // end of void  MultiLayerMieApplied::AddTargetLayer(...)  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::SetTargetPEC(double radius) {
    MarkUncalculated();
    if (target_width_.size() != 0 || target_index_.size() != 0)
      throw std::invalid_argument("Error! Define PEC target radius before any other layers!");
    // Add layer of any index...
    AddTargetLayer(radius, std::complex<double>(0.0, 0.0));
    // ... and mark it as PEC
    SetPECLayer(0);
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::SetCoatingIndex(std::vector<std::complex<double> > index) {
    MarkUncalculated();
    coating_index_.clear();
    for (auto value : index) coating_index_.push_back(value);
  }  // end of void MultiLayerMieApplied::SetCoatingIndex(std::vector<complex> index);  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::SetCoatingWidth(std::vector<double> width) {
    MarkUncalculated();
    coating_width_.clear();
    for (auto w : width)
      if (w <= 0)
        throw std::invalid_argument("Coating width should be positive!");
      else coating_width_.push_back(w);
  }
  // end of void MultiLayerMieApplied::SetCoatingWidth(...);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::SetWidthSP(const std::vector<double>& size_parameter) {
    MarkUncalculated();
    size_param_.clear();
    double prev_size_parameter = 0.0;
    for (auto layer_size_parameter : size_parameter) {
      if (layer_size_parameter <= 0.0)
        throw std::invalid_argument("Size parameter should be positive!");
      if (prev_size_parameter > layer_size_parameter) 
        throw std::invalid_argument
          ("Size parameter for next layer should be larger than the previous one!");
      prev_size_parameter = layer_size_parameter;
      size_param_.push_back(layer_size_parameter);
    }
  }
  // end of void MultiLayerMieApplied::SetWidthSP(...);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::SetIndexSP(const std::vector< std::complex<double> >& index) {
    MarkUncalculated();
    //refractive_index_.clear();
    refractive_index_ = index;
    // for (auto value : index) refractive_index_.push_back(value);
  }  // end of void MultiLayerMieApplied::SetIndexSP(...);  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::SetFieldPointsSP(const std::vector< std::vector<double> >& coords_sp) {
    if (coords_sp.size() != 3)
      throw std::invalid_argument("Error! Wrong dimension of field monitor points!");
    if (coords_sp[0].size() != coords_sp[1].size() || coords_sp[0].size() != coords_sp[2].size())
      throw std::invalid_argument("Error! Missing coordinates for field monitor points!");
    coords_sp_ = coords_sp;
    // for (int i = 0; i < coords_sp_[0].size(); ++i) {
    //   printf("%g, %g, %g\n", coords_sp_[0][i], coords_sp_[1][i], coords_sp_[2][i]);
    // }
  }  // end of void MultiLayerMieApplied::SetFieldPointsSP(...)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::GenerateSizeParameter() {
    MarkUncalculated();
    size_param_.clear();
    double radius = 0.0;
    for (auto width : target_width_) {
      radius += width;
      size_param_.push_back(2*PI_*radius/wavelength_);
    }
    for (auto width : coating_width_) {
      radius += width;
      size_param_.push_back(2*PI_*radius/wavelength_);
    }
    total_radius_ = radius;
  }  // end of void MultiLayerMieApplied::GenerateSizeParameter();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::GenerateIndex() {
    MarkUncalculated();
    refractive_index_.clear();
    for (auto index : target_index_) refractive_index_.push_back(index);
    for (auto index : coating_index_) refractive_index_.push_back(index);
  }  // end of void MultiLayerMieApplied::GenerateIndex();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double MultiLayerMieApplied::GetTotalRadius() {
    if (!isMieCalculated())  GenerateSizeParameter();
    return total_radius_;      
  }  // end of double MultiLayerMieApplied::GetTotalRadius();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector< std::vector<double> >
  MultiLayerMieApplied::GetSpectra(double from_WL, double to_WL, int samples) {
    if (!isMieCalculated())
      throw std::invalid_argument("You should run calculations before result request!");
    std::vector< std::vector<double> > spectra;
    double step_WL = (to_WL - from_WL)/static_cast<double>(samples);
    double wavelength_backup = wavelength_;
    long fails = 0;
    for (double WL = from_WL; WL < to_WL; WL += step_WL) {
      wavelength_ = WL;
      try {
        RunMieCalculation();
      } catch(const std::invalid_argument& ia) {
        fails++;
        continue;
      }
      //printf("%3.1f ",WL);
      spectra.push_back(std::vector<double>({wavelength_, GetQext(),
	      GetQsca(), GetQabs(), GetQbk()}));
    }  // end of for each WL in spectra
    printf("Spectrum has %li fails\n",fails);
    wavelength_ = wavelength_backup;
    return spectra;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::ClearTarget() {
    MarkUncalculated();
    target_width_.clear();
    target_index_.clear();
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::ClearCoating() {
    MarkUncalculated();
    coating_width_.clear();
    coating_index_.clear();
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::ClearLayers() {
    MarkUncalculated();
    ClearTarget();
    ClearCoating();
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::ClearAllDesign() {
    MarkUncalculated();
    ClearLayers();
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
  std::complex<double> MultiLayerMieApplied::calcD1confra(const int N, const std::complex<double> z) {
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
    CAK    = static_cast<std::complex<double> >(MM*KK)*ZINV; //debug -3 ZINV
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
      CAK = static_cast<std::complex<double> >(MM*KK)*ZINV; //    ** Eq. R25b //debug 5zinv
     //  //c ** Eq. R32    Ill-conditioned case -- stride two terms instead of one
     //  if (std::abs(CNUMER/CAK) >= EPS1 ||  std::abs(CDENOM/CAK) >= EPS1) {
     //         //c                       ** Eq. R34
     //         CNTN   = CAK*CNUMER + 1.0;
     //         CDTD   = CAK*CDENOM + 1.0;
     //         CONFRA = (CNTN/CDTD)*CONFRA; // ** Eq. R33
     //         MM  *= - 1;        KK  += 2;
     //         CAK = static_cast<std::complex<double> >(MM*KK)*ZINV; // ** Eq. R25b
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
       CONFRA = CAPT*CONFRA; // ** Eq. R26
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
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::ConvertToSP() {
    MarkUncalculated();
    if (target_width_.size() + coating_width_.size() == 0)
      return;  // Nothing to convert, we suppose that SP was set directly
    GenerateSizeParameter();
    GenerateIndex();
    if (size_param_.size() != refractive_index_.size())
      throw std::invalid_argument("Ivalid conversion of width to size parameter units!/n");
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::RunMieCalculation() {
    ConvertToSP();
    MultiLayerMie::RunMieCalculation(); 
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMieApplied::GetExpanCoeffs( std::vector< std::vector<std::complex<double> > >& aln, std::vector< std::vector<std::complex<double> > >& bln, std::vector< std::vector<std::complex<double> > >& cln, std::vector< std::vector<std::complex<double> > >& dln) {
    ConvertToSP();  // Case of call before running full Mie calculation.
    // Calculate scattering coefficients an_ and bn_
    calcScattCoeffs();
    // Calculate expansion coefficients aln_,  bln_, cln_, and dln_
    calcExpanCoeffs();
    aln = aln_;
    bln = bln_;
    cln = cln_;
    dln = dln_;
    
  }  // end of void MultiLayerMieApplied::GetExpanCoeffs( ...)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //

}  // end of namespace nmie
