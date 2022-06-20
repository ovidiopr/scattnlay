#ifndef SRC_NMIE_HPP_
#define SRC_NMIE_HPP_
//******************************************************************************
//    Copyright (C) 2009-2022  Ovidio Pena <ovidio@bytesfall.com>
//    Copyright (C) 2013-2022  Konstantin Ladutenko <kostyfisik@gmail.com>
//
//    This file is part of scattnlay
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    The only additional remark is that we expect that all publications
//    describing work using this software, or all commercial products
//    using it, cite at least one of the following references:
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by
//        a multilayered sphere," Computer Physics Communications,
//        vol. 180, Nov. 2009, pp. 2348-2354.
//    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie
//        calculation of electromagnetic near-field for a multilayered
//        sphere," Computer Physics Communications, vol. 214, May 2017,
//        pp. 225-230.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//******************************************************************************

#define VERSION "2.2"  // Compare with Makefile and setup.py
#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "nmie-precision.hpp"
//#ifdef MULTI_PRECISION
//#include <boost/math/constants/constants.hpp>
//#endif

namespace nmie {
//******************************************************************************
int ScattCoeffs(const unsigned int L,
                const int pl,
                const std::vector<double>& x,
                const std::vector<std::complex<double>>& m,
                const int nmax,
                std::vector<std::complex<double>>& an,
                std::vector<std::complex<double>>& bn);

int ExpanCoeffs(const unsigned int L,
                const int pl,
                const std::vector<double>& x,
                const std::vector<std::complex<double>>& m,
                const int nmax,
                std::vector<std::vector<std::complex<double>>>& an,
                std::vector<std::vector<std::complex<double>>>& bn,
                std::vector<std::vector<std::complex<double>>>& cn,
                std::vector<std::vector<std::complex<double>>>& dn);

//******************************************************************************
// helper functions
//******************************************************************************
template <typename FloatType>
double eval_delta(const unsigned int steps,
                  const double from_value,
                  const double to_value);

template <class T>
inline T pow2(const T value) {
  return value * value;
}

template <class T>
inline T cabs(const std::complex<T> value) {
  return sqrt_t(pow2(value.real()) + pow2(value.imag()));
}

template <class T>
inline T vabs(const std::vector<std::complex<T>> value) {
  return nmm::sqrt(pow2(value[0].real()) + pow2(value[1].real()) +
                   pow2(value[2].real()) + pow2(value[0].imag()) +
                   pow2(value[1].imag()) + pow2(value[2].imag()));
}

template <typename FloatType>
int newround(FloatType x) {
  return x >= 0 ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
  // return x >= 0 ? (x + 0.5).convert_to<int>():(x - 0.5).convert_to<int>();
}

template <typename T>
inline std::complex<T> my_exp(const std::complex<T>& x) {
  using std::exp;  // use ADL
  T const& r = exp(x.real());
  return std::polar(r, x.imag());
}

//******************************************************************************
// pl, nmax, mode_n, mode_type
int nMie(const unsigned int L,
         const int pl,
         std::vector<double>& x,
         std::vector<std::complex<double>>& m,
         const unsigned int nTheta,
         std::vector<double>& Theta,
         const int nmax,
         double* Qext,
         double* Qsca,
         double* Qabs,
         double* Qbk,
         double* Qpr,
         double* g,
         double* Albedo,
         std::vector<std::complex<double>>& S1,
         std::vector<std::complex<double>>& S2,
         int mode_n,
         int mode_type);

//******************************************************************************
// pl and nmax
int nMie(const unsigned int L,
         const int pl,
         std::vector<double>& x,
         std::vector<std::complex<double>>& m,
         const unsigned int nTheta,
         std::vector<double>& Theta,
         const int nmax,
         double* Qext,
         double* Qsca,
         double* Qabs,
         double* Qbk,
         double* Qpr,
         double* g,
         double* Albedo,
         std::vector<std::complex<double>>& S1,
         std::vector<std::complex<double>>& S2);

//******************************************************************************
// no pl and nmax
int nMie(const unsigned int L,
         std::vector<double>& x,
         std::vector<std::complex<double>>& m,
         const unsigned int nTheta,
         std::vector<double>& Theta,
         double* Qext,
         double* Qsca,
         double* Qabs,
         double* Qbk,
         double* Qpr,
         double* g,
         double* Albedo,
         std::vector<std::complex<double>>& S1,
         std::vector<std::complex<double>>& S2);

//******************************************************************************
// pl
int nMie(const unsigned int L,
         const int pl,
         std::vector<double>& x,
         std::vector<std::complex<double>>& m,
         const unsigned int nTheta,
         std::vector<double>& Theta,
         double* Qext,
         double* Qsca,
         double* Qabs,
         double* Qbk,
         double* Qpr,
         double* g,
         double* Albedo,
         std::vector<std::complex<double>>& S1,
         std::vector<std::complex<double>>& S2);

//******************************************************************************
// nmax
int nMie(const unsigned int L,
         std::vector<double>& x,
         std::vector<std::complex<double>>& m,
         const unsigned int nTheta,
         std::vector<double>& Theta,
         const int nmax,
         double* Qext,
         double* Qsca,
         double* Qabs,
         double* Qbk,
         double* Qpr,
         double* g,
         double* Albedo,
         std::vector<std::complex<double>>& S1,
         std::vector<std::complex<double>>& S2);

//******************************************************************************
int nField(const unsigned int L,
           const int pl,
           const std::vector<double>& x,
           const std::vector<std::complex<double>>& m,
           const int nmax,
           const int mode_n,
           const int mode_type,
           const unsigned int ncoord,
           const std::vector<double>& Xp,
           const std::vector<double>& Yp,
           const std::vector<double>& Zp,
           std::vector<std::vector<std::complex<double>>>& E,
           std::vector<std::vector<std::complex<double>>>& H);

//******************************************************************************
// constants for per mode evaluation
//******************************************************************************
enum Modes { kAll = -1, kElectric = 0, kMagnetic = 1 };
enum Planes { kEk = 0, kHk = 1, kEH = 2 };

//******************************************************************************
const FloatType PI_ =
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949129833673362440656643086021394946395224737190702179860943702770539217176293176752384674818467669405132000568127145263560827785771342757789609173637178721468440901224953430146549585371050792279689258923542019956112129021960864034418159813629774771309960518707211349999998372978049951059731732816096318595024459455346908302642522308253344685035261931188171010003137838752886587533208381420617177669147303598253490428755468731159562863882353787593751957781857780532171226806613001927876611195909216420198938095257201065485863278865936153381827968230301952035301852968995773622599413891249721775283479131515574857242454150695950829533116861727855889075098381754637464939319255060400927701671139009848824012858361603563707660104710181942955596198946767837449448255379774726847104047534646208046684259069491293313677028989152104752162056966024058038150193511253382430035587640247496473263914199272604269922796782354781636009341721641219924586315030286182974555706749838505494588586926995690927210797509302955321165344987202755960236480665499119881834797753566369807426542527862551818417574672890977772793800081647060016145249192173217214772350141441973568548161361157352552133475741849468438523323907394143334547762416862518983569485562099219222184272550254256887671790494601653466804988627232791786085784383827967976681454100953883786360950680064225125205117392984896084128488626945604241965285022210661186306744278622039194945047123713786960956364371917287467764657573962413890865832645995813390478027590099465764078951269468398352595709825822620522489407726719478268482601476990902640136394437455305068203496252451749399651431429809190659250937221696461515709858387410597885959772975498930161753928468138268683868942774155991855925245953959431049972524680845987273644695848653836736222626099124608051243884390451244136549762780797715691435997700129616089441694868555848406353422072225828488648158456028506016842739452267467678895252138522549954666727823986456596116354886230577456498035593634568174324112515076069479451096596094025228879710893145669136867228748940560101503308;
// light speed [m/s]
const double cc_ = 2.99792458e8;
// assume non-magnetic (MU=MU0=const) [N/A^2]
const FloatType mu_ = 4.0 * PI_ * 1.0e-7;

//******************************************************************************
//******************************************************************************
template <typename FloatType = double>
class MultiLayerMie {
 public:
#ifdef MULTI_PRECISION
  const FloatType convergence_threshold_ = std::pow(10, -MULTI_PRECISION / 2);
  //    const FloatType convergence_threshold_ = 1e-50;
  //    const FloatType nearfield_convergence_threshold_ = std::pow(10,
  //    -MULTI_PRECISION/2);

  // For near-field evaluation we use Le Ru cutoff which is valid only for
  // double precision, so convergence threshold is the same
  const FloatType nearfield_convergence_threshold_ = 1e-14;
#else
  const double convergence_threshold_ = 1e-25;
  const double nearfield_convergence_threshold_ = 1e-14;
#endif
  void RunMieCalculation();

  void RunFieldCalculation(bool isMarkUnconverged = true);

  void RunFieldCalculationPolar(
      const int outer_arc_points = 1,
      const int radius_points = 1,
      const double from_Rho = 0,
      const double to_Rho = static_cast<double>(1.),
      const double from_Theta = 0,
      const double to_Theta = static_cast<double>(3.14159265358979323),
      const double from_Phi = 0,
      const double to_Phi = static_cast<double>(3.14159265358979323),
      const bool isMarkUnconverged = true,
      int nmax_in = -1);

  void RunFieldCalculationCartesian(const int first_side_points = 2,
                                    const int second_side_points = 2,
                                    const double relative_side_length = 2,
                                    const int plane_selected = Planes::kEk,
                                    const double at_x = 0,
                                    const double at_y = 0,
                                    const double at_z = 0,
                                    const bool isMarkUnconverged = true,
                                    const int nmax_in = -1);

  void calcScattCoeffs();
  void calcExpanCoeffs();
  //****************************************************************************
  // Return calculation results
  //****************************************************************************
  template <typename outputType = FloatType>
  outputType GetQext();

  template <typename outputType = FloatType>
  outputType GetQsca();

  template <typename outputType = FloatType>
  outputType GetQabs();

  template <typename outputType = FloatType>
  outputType GetQbk();

  template <typename outputType = FloatType>
  outputType GetQpr();

  template <typename outputType = FloatType>
  outputType GetAsymmetryFactor();

  template <typename outputType = FloatType>
  outputType GetAlbedo();

  std::vector<std::complex<FloatType>> GetS1();
  std::vector<std::complex<FloatType>> GetS2();

  std::vector<std::complex<FloatType>> GetAn() {
    return an_;
  };
  std::vector<std::complex<FloatType>> GetBn() {
    return bn_;
  };

  std::vector<std::vector<std::complex<FloatType>>> GetLayerAn() {
    return aln_;
  };
  std::vector<std::vector<std::complex<FloatType>>> GetLayerBn() {
    return bln_;
  };
  std::vector<std::vector<std::complex<FloatType>>> GetLayerCn() {
    return cln_;
  };
  std::vector<std::vector<std::complex<FloatType>>> GetLayerDn() {
    return dln_;
  };

  //****************************************************************************
  // Problem definition
  // Modify size of all layers
  //****************************************************************************
  void SetLayersSize(const std::vector<FloatType>& layer_size);
  // Modify refractive index of all layers
  void SetLayersIndex(const std::vector<std::complex<FloatType>>& index);

  template <typename evalType = FloatType>
  void GetIndexAtRadius(const evalType Rho,
                        std::complex<evalType>& ml,
                        unsigned int& l);

  template <typename evalType = FloatType>
  void GetIndexAtRadius(const evalType Rho, std::complex<evalType>& ml);

  // Modify scattering (theta) angles
  void SetAngles(const std::vector<FloatType>& angles);

  // Modify coordinates for field calculation
  void SetFieldCoords(const std::vector<std::vector<FloatType>>& coords);

  // Modify index of PEC layer
  void SetPECLayer(int layer_position = 0);

  // Modify the mode taking into account for evaluation of output variables
  void SetModeNmaxAndType(int mode_n, int mode_type) {
    mode_n_ = mode_n;
    mode_type_ = mode_type;
  };

  // Set a fixed value for the maximum number of terms
  void SetMaxTerms(int nmax);

  // Get maximum number of terms
  int GetMaxTerms() {
    return nmax_;
  };

  bool isMieCalculated() {
    return isMieCalculated_;
  };

  // Clear layer information
  void ClearLayers();

  void MarkUncalculated();

  // Read parameters
  // Get total size parameter of particle
  FloatType GetSizeParameter();
  // Returns size of all layers
  std::vector<FloatType> GetLayersSize() {
    return size_param_;
  };
  // Returns refractive index of all layers
  std::vector<std::complex<FloatType>> GetLayersIndex() {
    return refractive_index_;
  };
  // Returns scattering (theta) angles
  std::vector<FloatType> GetAngles() {
    return theta_;
  };
  // Returns coordinates used for field calculation
  std::vector<std::vector<FloatType>> GetFieldCoords() {
    return coords_;
  };
  // Returns index of PEC layer
  int GetPECLayer() {
    return PEC_layer_position_;
  };

  std::vector<std::vector<std::complex<FloatType>>> GetFieldE() {
    return E_;
  };  // {X[], Y[], Z[]}
  std::vector<std::vector<std::complex<FloatType>>> GetFieldH() {
    return H_;
  };

  std::vector<FloatType> GetFieldEabs() {
    return Eabs_;
  };  // {X[], Y[], Z[]}
  std::vector<FloatType> GetFieldHabs() {
    return Habs_;
  };
  bool GetFieldConvergence();

  // Get fields in spherical coordinates.
  std::vector<std::vector<std::complex<FloatType>>> GetFieldEs() {
    return Es_;
  };  // {rho[], theta[], phi[]}
  std::vector<std::vector<std::complex<FloatType>>> GetFieldHs() {
    return Hs_;
  };

 protected:
  // Size parameter for all layers
  std::vector<FloatType> size_param_;
  // Refractive index for all layers
  std::vector<std::complex<FloatType>> refractive_index_;
  // Scattering coefficients
  std::vector<std::complex<FloatType>> an_, bn_;
  std::vector<std::vector<std::complex<FloatType>>> aln_, bln_, cln_, dln_;
  // Points for field evaluation
  std::vector<std::vector<FloatType>> coords_;
  std::vector<std::vector<FloatType>> coords_polar_;

 private:
  unsigned int calcNstop(FloatType xL = -1);
  unsigned int calcNmax(FloatType xL = -1);

  std::complex<FloatType> calc_an(int n,
                                  FloatType XL,
                                  std::complex<FloatType> Ha,
                                  std::complex<FloatType> mL,
                                  std::complex<FloatType> PsiXL,
                                  std::complex<FloatType> ZetaXL,
                                  std::complex<FloatType> PsiXLM1,
                                  std::complex<FloatType> ZetaXLM1);

  std::complex<FloatType> calc_bn(int n,
                                  FloatType XL,
                                  std::complex<FloatType> Hb,
                                  std::complex<FloatType> mL,
                                  std::complex<FloatType> PsiXL,
                                  std::complex<FloatType> ZetaXL,
                                  std::complex<FloatType> PsiXLM1,
                                  std::complex<FloatType> ZetaXLM1);

  std::complex<FloatType> calc_S1(int n,
                                  std::complex<FloatType> an,
                                  std::complex<FloatType> bn,
                                  FloatType Pi,
                                  FloatType Tau);

  std::complex<FloatType> calc_S2(int n,
                                  std::complex<FloatType> an,
                                  std::complex<FloatType> bn,
                                  FloatType Pi,
                                  FloatType Tau);

  void calcD1D3(std::complex<FloatType> z,
                std::vector<std::complex<FloatType>>& D1,
                std::vector<std::complex<FloatType>>& D3);

  void calcPsiZeta(std::complex<FloatType> x,
                   std::vector<std::complex<FloatType>>& Psi,
                   std::vector<std::complex<FloatType>>& Zeta);
  void calcPiTau(const FloatType& costheta,
                 std::vector<FloatType>& Pi,
                 std::vector<FloatType>& Tau);

  template <typename evalType = FloatType>
  void calcSpherHarm(const std::complex<evalType> Rho,
                     const evalType Theta,
                     const evalType Phi,
                     const std::complex<evalType>& rn,
                     const std::complex<evalType>& Dn,
                     const evalType& Pi,
                     const evalType& Tau,
                     const evalType& n,
                     std::vector<std::complex<evalType>>& Mo1n,
                     std::vector<std::complex<evalType>>& Me1n,
                     std::vector<std::complex<evalType>>& No1n,
                     std::vector<std::complex<evalType>>& Ne1n);

  template <typename evalType = FloatType>
  void calcFieldByComponents(const evalType Rho,
                             const evalType Theta,
                             const evalType Phi,
                             const std::vector<std::complex<evalType>>& Psi,
                             const std::vector<std::complex<evalType>>& D1n,
                             const std::vector<std::complex<evalType>>& Zeta,
                             const std::vector<std::complex<evalType>>& D3n,
                             const std::vector<evalType>& Pi,
                             const std::vector<evalType>& Tau,
                             std::vector<std::complex<evalType>>& E,
                             std::vector<std::complex<evalType>>& H,
                             std::vector<bool>& isConvergedE,
                             std::vector<bool>& isConvergedH,
                             bool isMarkUnconverged);

  bool isExpCoeffsCalc_ = false;
  bool isScaCoeffsCalc_ = false;
  bool isMieCalculated_ = false;
  std::vector<bool> isConvergedE_ = {false, false, false};
  std::vector<bool> isConvergedH_ = {false, false, false};

  // Scattering angles for scattering pattern in radians
  std::vector<FloatType> theta_;

  // Should be -1 if there is no PEC.
  int PEC_layer_position_ = -1;

  int mode_n_ = Modes::kAll;
  int mode_type_ = Modes::kAll;

  // with calcNmax(int first_layer);
  int nmax_ = -1;
  int nmax_preset_ = -1;
  int available_maximal_nmax_ = -1;

  // Store result
  FloatType Qsca_ = 0.0, Qext_ = 0.0, Qabs_ = 0.0, Qbk_ = 0.0, Qpr_ = 0.0,
            asymmetry_factor_ = 0.0, albedo_ = 0.0;

  // {X[], Y[], Z[]}
  std::vector<std::vector<std::complex<FloatType>>> E_, H_;
  std::vector<std::vector<std::complex<FloatType>>> Es_, Hs_;
  std::vector<FloatType> Eabs_, Habs_;

  std::vector<std::complex<FloatType>> S1_, S2_;
  void calcMieSeriesNeededToConverge(const FloatType Rho, int nmax_in = -1);

  void calcPiTauAllTheta(const double from_Theta,
                         const double to_Theta,
                         std::vector<std::vector<FloatType>>& Pi,
                         std::vector<std::vector<FloatType>>& Tau);

  void calcRadialOnlyDependantFunctions(
      const double from_Rho,
      const double to_Rho,
      std::vector<std::vector<std::complex<FloatType>>>& Psi,
      std::vector<std::vector<std::complex<FloatType>>>& D1n,
      std::vector<std::vector<std::complex<FloatType>>>& Zeta,
      std::vector<std::vector<std::complex<FloatType>>>& D3n,
      int nmax_in = -1);

  void convertFieldsFromSphericalToCartesian();

  void UpdateConvergenceStatus(std::vector<bool> isConvergedE,
                               std::vector<bool> isConvergedH);
};  // end of class MultiLayerMie

//******************************************************************************
//******************************************************************************
template <typename FloatType = double>
class MesoMie {
 public:
  std::vector<std::complex<FloatType>> an_, bn_;
  FloatType x_;
  std::complex<FloatType> m_;
  std::vector<std::complex<FloatType>> GetAn() { return an_; };
  std::vector<std::complex<FloatType>> GetBn() { return bn_; };

  FloatType Qsca_ = 0.0, Qext_ = 0.0;

  template <typename outputType = FloatType>
  outputType GetQsca() {
    return static_cast<outputType>(Qsca_);
  }

  template <typename outputType = FloatType>
  outputType GetQext() {
    return static_cast<outputType>(Qext_);
  }

  void calc_ab(FloatType R,
               FloatType xd,
               std::complex<FloatType> xm,
               std::complex<FloatType> eps_d,
               std::complex<FloatType> eps_m,
               std::complex<FloatType> d_parallel,
               std::complex<FloatType> d_perp);

  void calc_Q();
  // template <typename outputType = FloatType>
  // outputType GetQext();

  // template <typename outputType = FloatType>
  // outputType GetQsca();
};  // end of class MesoMie

}  // end of namespace nmie
#endif  // SRC_NMIE_HPP_
