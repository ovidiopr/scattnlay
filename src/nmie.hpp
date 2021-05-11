#ifndef SRC_NMIE_HPP_
#define SRC_NMIE_HPP_
//**********************************************************************************//
//    Copyright (C) 2009-2018  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2018  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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
//    using it, cite at least one of the following references:                      //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie              //
//        calculation of electromagnetic near-field for a multilayered              //
//        sphere," Computer Physics Communications, vol. 214, May 2017,             //
//        pp. 225-230.                                                              //
//                                                                                  //
//    You should have received a copy of the GNU General Public License             //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//**********************************************************************************//

#define VERSION "2.2"  //Compare with Makefile and setup.py
#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "nmie-precision.hpp"
#ifdef MULTI_PRECISION
#include <boost/math/constants/constants.hpp>
#endif
namespace py = pybind11;

template <typename T>
std::vector<T> Py2Vector(const py::array_t<T> &py_x) {
  std::vector<T> c_x(py_x.size());
  std::memcpy(c_x.data(), py_x.data(), py_x.size()*sizeof(T));
  return c_x;
}

// https://github.com/pybind/pybind11/issues/1042#issuecomment-508582847
//template <typename Sequence>
//inline py::array_t<typename Sequence::value_type> Vector2Py(Sequence&& seq) {
//  // Move entire object to heap (Ensure is moveable!). Memory handled via Python capsule
//  Sequence* seq_ptr = new Sequence(std::move(seq));
//  auto capsule = py::capsule(seq_ptr, [](void* p) { delete reinterpret_cast<Sequence*>(p); });
//  return py::array(seq_ptr->size(),  // shape of array
//                   seq_ptr->data(),  // c-style contiguous strides for Sequence
//                   capsule           // numpy array references this parent
//  );
//}
template <typename outputType>
inline py::array_t<outputType> Vector2Py(const std::vector<outputType>& seq) {
  return py::array(seq.size(), seq.data());
}

template <typename inputType=double, typename outputType=double>
py::array_t< std::complex<outputType>> VectorComplex2Py(const std::vector<std::complex<inputType> > &cf_x) {
  auto c_x = nmie::ConvertComplexVector<outputType, inputType>(cf_x);
  auto py_x = py::array_t< std::complex<outputType>>(c_x.size());
  auto py_x_buffer = py_x.request();
  auto *py_x_ptr = (std::complex<outputType> *) py_x_buffer.ptr;
  std::memcpy(py_x_ptr, c_x.data(), c_x.size()*sizeof(std::complex<outputType>));
  return py_x;
}

// https://stackoverflow.com/questions/17294629/merging-flattening-sub-vectors-into-a-single-vector-c-converting-2d-to-1d
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &v) {
  std::size_t total_size = 0;
  for (const auto &sub : v)
    total_size += sub.size(); // I wish there was a transform_accumulate
  std::vector<T> result;
  result.reserve(total_size);
  for (const auto &sub : v)
    result.insert(result.end(), sub.begin(), sub.end());
  return result;
}


template <typename T>
py::array Vector2DComplex2Py(const std::vector<std::vector<T > > &x) {
  size_t dim1 = x.size();
  size_t dim2 = x[0].size();
  auto result = flatten(x);
  // https://github.com/tdegeus/pybind11_examples/blob/master/04_numpy-2D_cpp-vector/example.cpp
  size_t              ndim    = 2;
  std::vector<size_t> shape   = {dim1, dim2};
  std::vector<size_t> strides = {sizeof(T)*dim2, sizeof(T)};

  // return 2-D NumPy array
  return py::array(py::buffer_info(
      result.data(),                       /* data as contiguous array  */
      sizeof(T),                           /* size of one scalar        */
      py::format_descriptor<T>::format(),  /* data type                 */
      ndim,                                /* number of dimensions      */
      shape,                               /* shape of the matrix       */
      strides                              /* strides for each axis     */
  ));
}

namespace nmie {
  int ScattCoeffs(const unsigned int L, const int pl,
                  const std::vector<double> &x, const std::vector<std::complex<double> > &m,
                  const int nmax,
                  std::vector<std::complex<double> > &an,
                  std::vector<std::complex<double> > &bn);

  int ExpanCoeffs(const unsigned int L, const int pl,
                  const std::vector<double> &x, const std::vector<std::complex<double> > &m,
                  const int nmax,
                  std::vector<std::vector<std::complex<double> > > &an,
                  std::vector<std::vector<std::complex<double> > > &bn,
                  std::vector<std::vector<std::complex<double> > > &cn,
                  std::vector<std::vector<std::complex<double> > > &dn);

//helper functions
template <typename FloatType>
double eval_delta(const unsigned int steps, const double from_value, const double to_value);


template<class T> inline T pow2(const T value) {return value*value;}

template<class T> inline T cabs(const std::complex<T> value)
{return nmm::sqrt(pow2(value.real()) + pow2(value.imag()));}

template<class T> inline T vabs(const std::vector<std::complex<T>> value)
{return nmm::sqrt(
    pow2(value[0].real()) + pow2(value[1].real()) + pow2(value[2].real())
    +pow2(value[0].imag()) + pow2(value[1].imag()) + pow2(value[2].imag()));}

template <typename FloatType>
int newround(FloatType x) {
  return x >= 0 ? static_cast<int>(x + 0.5):static_cast<int>(x - 0.5);
  //return x >= 0 ? (x + 0.5).convert_to<int>():(x - 0.5).convert_to<int>();
}
template<typename T>
inline std::complex<T> my_exp(const std::complex<T> &x) {
  using std::exp; // use ADL
  T const &r = exp(x.real());
  return std::polar(r, x.imag());
}

// pl, nmax, mode_n, mode_type
    int nMie(const unsigned int L,
           const int pl,
           std::vector<double> &x, std::vector<std::complex<double> > &m,
           const unsigned int nTheta, std::vector<double> &Theta,
           const int nmax,
           double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr,
           double *g, double *Albedo,
           std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2,
           int mode_n, int mode_type);
  // pl and nmax
    int nMie(const unsigned int L,
           const int pl,
           std::vector<double> &x, std::vector<std::complex<double> > &m,
           const unsigned int nTheta, std::vector<double> &Theta,
           const int nmax,
           double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr,
           double *g, double *Albedo,
           std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2);
  // no pl and nmax
  int nMie(const unsigned int L,
           std::vector<double> &x, std::vector<std::complex<double> > &m,
           const unsigned int nTheta, std::vector<double> &Theta,
           double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr,
           double *g, double *Albedo,
           std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2);
  // pl
  int nMie(const unsigned int L,
           const int pl,
           std::vector<double> &x, std::vector<std::complex<double> > &m,
           const unsigned int nTheta, std::vector<double> &Theta,
           double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr,
           double *g, double *Albedo,
           std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2);
  // nmax
  int nMie(const unsigned int L,
           std::vector<double> &x, std::vector<std::complex<double> > &m,
           const unsigned int nTheta, std::vector<double> &Theta,
           const int nmax,
           double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr,
           double *g, double *Albedo,
           std::vector<std::complex<double> > &S1, std::vector<std::complex<double> > &S2);
  int nField(const unsigned int L, const int pl,
             const std::vector<double> &x, const std::vector<std::complex<double> > &m, const int nmax,
             const int mode_n, const int mode_type, const unsigned int ncoord,
             const std::vector<double> &Xp, const std::vector<double> &Yp, const std::vector<double> &Zp,
             std::vector<std::vector<std::complex<double> > > &E, std::vector<std::vector<std::complex<double> > > &H);

  // constants for per mode evaluation
  enum Modes {kAll = -1, kElectric = 0, kMagnetic = 1};

  template <typename FloatType = double>
  class MultiLayerMie {
   public:
    const FloatType PI_=3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949129833673362440656643086021394946395224737190702179860943702770539217176293176752384674818467669405132000568127145263560827785771342757789609173637178721468440901224953430146549585371050792279689258923542019956112129021960864034418159813629774771309960518707211349999998372978049951059731732816096318595024459455346908302642522308253344685035261931188171010003137838752886587533208381420617177669147303598253490428755468731159562863882353787593751957781857780532171226806613001927876611195909216420198938095257201065485863278865936153381827968230301952035301852968995773622599413891249721775283479131515574857242454150695950829533116861727855889075098381754637464939319255060400927701671139009848824012858361603563707660104710181942955596198946767837449448255379774726847104047534646208046684259069491293313677028989152104752162056966024058038150193511253382430035587640247496473263914199272604269922796782354781636009341721641219924586315030286182974555706749838505494588586926995690927210797509302955321165344987202755960236480665499119881834797753566369807426542527862551818417574672890977772793800081647060016145249192173217214772350141441973568548161361157352552133475741849468438523323907394143334547762416862518983569485562099219222184272550254256887671790494601653466804988627232791786085784383827967976681454100953883786360950680064225125205117392984896084128488626945604241965285022210661186306744278622039194945047123713786960956364371917287467764657573962413890865832645995813390478027590099465764078951269468398352595709825822620522489407726719478268482601476990902640136394437455305068203496252451749399651431429809190659250937221696461515709858387410597885959772975498930161753928468138268683868942774155991855925245953959431049972524680845987273644695848653836736222626099124608051243884390451244136549762780797715691435997700129616089441694868555848406353422072225828488648158456028506016842739452267467678895252138522549954666727823986456596116354886230577456498035593634568174324112515076069479451096596094025228879710893145669136867228748940560101503308;
    // light speed [m s-1]
    const double cc_ = 2.99792458e8;
    // assume non-magnetic (MU=MU0=const) [N A-2]
    const FloatType mu_ = 4.0*PI_*1.0e-7;
    // Run calculation
    void RunMieCalculation();
    void RunFieldCalculation();
    void RunFieldCalculationPolar(const int outer_arc_points = 1,
                                  const int radius_points=1,
                                  const double from_Rho=0, const double to_Rho=static_cast<double>(1.),
                                  const double from_Theta=0, const double to_Theta=static_cast<double>(3.14159265358979323),
                                  const double from_Phi=0, const double to_Phi=static_cast<double>(3.14159265358979323),
                                  const bool isIgnoreAvailableNmax = false);

    void calcScattCoeffs();
    void calcExpanCoeffs();

    // Return calculation results
    template <typename outputType = FloatType> outputType GetQext();
    template <typename outputType = FloatType> outputType  GetQsca();
    template <typename outputType = FloatType> outputType  GetQabs();
    template <typename outputType = FloatType> outputType  GetQbk();
    template <typename outputType = FloatType> outputType  GetQpr();
    template <typename outputType = FloatType> outputType  GetAsymmetryFactor();
    template <typename outputType = FloatType> outputType  GetAlbedo();
    std::vector<std::complex<FloatType> > GetS1();
    std::vector<std::complex<FloatType> > GetS2();
    template <typename outputType> py::array_t< std::complex<outputType>>  GetS1();
    template <typename outputType> py::array_t< std::complex<outputType>>  GetS2();

    std::vector<std::complex<FloatType> > GetAn(){return an_;};
    std::vector<std::complex<FloatType> > GetBn(){return bn_;};
    template <typename outputType> py::array_t< std::complex<outputType>>  GetAn();
    template <typename outputType> py::array_t< std::complex<outputType>>  GetBn();

    std::vector< std::vector<std::complex<FloatType> > > GetLayerAn(){return aln_;};
    std::vector< std::vector<std::complex<FloatType> > > GetLayerBn(){return bln_;};
    std::vector< std::vector<std::complex<FloatType> > > GetLayerCn(){return cln_;};
    std::vector< std::vector<std::complex<FloatType> > > GetLayerDn(){return dln_;};
    template <typename outputType> py::array GetLayerAn();
    template <typename outputType> py::array GetLayerBn();
    template <typename outputType> py::array GetLayerCn();
    template <typename outputType> py::array GetLayerDn();

    // Problem definition
    // Modify size of all layers
    void SetLayersSize(const std::vector<FloatType> &layer_size);
    template <typename inputType>
    void SetLayersSize(const py::array_t<inputType, py::array::c_style | py::array::forcecast> &py_layer_size);
    // Modify refractive index of all layers
    void SetLayersIndex(const std::vector< std::complex<FloatType> > &index);
    template <typename inputType>
    void SetLayersIndex(const py::array_t<std::complex<inputType>, py::array::c_style | py::array::forcecast> &py_index);

    template <typename evalType=FloatType> void GetIndexAtRadius(const evalType Rho, std::complex<evalType> &ml, unsigned int &l);
    template <typename evalType=FloatType> void GetIndexAtRadius(const evalType Rho, std::complex<evalType> &ml);
    // Modify scattering (theta) py_angles
    void SetAngles(const std::vector<FloatType> &py_angles);
    template <typename inputType>
    void SetAngles(const py::array_t<inputType, py::array::c_style | py::array::forcecast> &py_angles);
    // Modify coordinates for field calculation
    void SetFieldCoords(const std::vector< std::vector<FloatType> > &coords);
    void SetFieldCoords(const py::array_t<double, py::array::c_style | py::array::forcecast> &py_Xp,
                        const py::array_t<double, py::array::c_style | py::array::forcecast> &py_Yp,
                        const py::array_t<double, py::array::c_style | py::array::forcecast> &py_Zp);
    // Modify index of PEC layer
    void SetPECLayer(int layer_position = 0);
    // Modify the mode taking into account for evaluation of output variables
    void SetModeNmaxAndType(int mode_n, int mode_type){mode_n_ = mode_n; mode_type_ = mode_type;};

    // Set a fixed value for the maximun number of terms
    void SetMaxTerms(int nmax);
    // Get maximum number of terms
    int GetMaxTerms() {return nmax_;};

    bool isMieCalculated(){return isMieCalculated_;};
    // Clear layer information
    void ClearLayers();
    void MarkUncalculated();

    // Read parameters
    // Get total size parameter of particle
    FloatType GetSizeParameter();
    // Returns size of all layers
    std::vector<FloatType> GetLayersSize(){return size_param_;};
    // Returns refractive index of all layers
    std::vector<std::complex<FloatType> > GetLayersIndex(){return refractive_index_;};
    // Returns scattering (theta) angles
    std::vector<FloatType> GetAngles(){return theta_;};
    // Returns coordinates used for field calculation
    std::vector<std::vector<FloatType> > GetFieldCoords(){return coords_;};
    // Returns index of PEC layer
    int GetPECLayer(){return PEC_layer_position_;};

    std::vector<std::vector< std::complex<FloatType> > > GetFieldE(){return E_;};   // {X[], Y[], Z[]}
    std::vector<std::vector< std::complex<FloatType> > > GetFieldH(){return H_;};
    template <typename outputType> py::array GetFieldE();
    template <typename outputType> py::array GetFieldH();

    std::vector< FloatType> GetFieldEabs(){return Eabs_;};   // {X[], Y[], Z[]}
    std::vector< FloatType> GetFieldHabs(){return Habs_;};
//    template <typename outputType> py::array_t<outputType>  GetFieldEabs();
//    template <typename outputType> py::array_t<outputType>  GetFieldHabs();

    // Get fields in spherical coordinates.
    std::vector<std::vector< std::complex<FloatType> > > GetFieldEs(){return Es_;};   // {rho[], teha[], phi[]}
    std::vector<std::vector< std::complex<FloatType> > > GetFieldHs(){return Hs_;};

  protected:
    // Size parameter for all layers
    std::vector<FloatType> size_param_;
    // Refractive index for all layers
    std::vector< std::complex<FloatType> > refractive_index_;
    // Scattering coefficients
    std::vector<std::complex<FloatType> > an_, bn_;
    std::vector< std::vector<std::complex<FloatType> > > aln_, bln_, cln_, dln_;
    // Points for field evaluation
    std::vector< std::vector<FloatType> > coords_;
    std::vector< std::vector<FloatType> > coords_polar_;

  private:
    unsigned int calcNstop(FloatType xL = -1);
    unsigned int calcNmax(FloatType xL = -1);

    std::complex<FloatType> calc_an(int n, FloatType XL, std::complex<FloatType> Ha, std::complex<FloatType> mL,
                                 std::complex<FloatType> PsiXL, std::complex<FloatType> ZetaXL,
                                 std::complex<FloatType> PsiXLM1, std::complex<FloatType> ZetaXLM1);
    std::complex<FloatType> calc_bn(int n, FloatType XL, std::complex<FloatType> Hb, std::complex<FloatType> mL,
                                 std::complex<FloatType> PsiXL, std::complex<FloatType> ZetaXL,
                                 std::complex<FloatType> PsiXLM1, std::complex<FloatType> ZetaXLM1);
    std::complex<FloatType> calc_S1(int n, std::complex<FloatType> an, std::complex<FloatType> bn,
                                 FloatType Pi, FloatType Tau);
    std::complex<FloatType> calc_S2(int n, std::complex<FloatType> an, std::complex<FloatType> bn,
                                 FloatType Pi, FloatType Tau);
    void calcD1D3(std::complex<FloatType> z,
                  std::vector<std::complex<FloatType> > &D1,
                  std::vector<std::complex<FloatType> > &D3);
    void calcPsiZeta(std::complex<FloatType> x,
                     std::vector<std::complex<FloatType> > &Psi,
                     std::vector<std::complex<FloatType> > &Zeta);
    void calcPiTau(const FloatType &costheta,
                   std::vector<FloatType> &Pi, std::vector<FloatType> &Tau);
    template <typename evalType=FloatType>
    void calcSpherHarm(const std::complex<evalType> Rho, const evalType Theta, const evalType Phi,
                       const std::complex<evalType> &rn, const std::complex<evalType> &Dn,
                       const evalType &Pi, const evalType &Tau, const evalType &n,
                       std::vector<std::complex<evalType> > &Mo1n, std::vector<std::complex<evalType> > &Me1n,
                       std::vector<std::complex<evalType> > &No1n, std::vector<std::complex<evalType> > &Ne1n);

    template <typename evalType=FloatType>
    void calcFieldByComponents(const evalType Rho, const evalType Theta, const evalType Phi,
                               const std::vector<std::complex<evalType> > &Psi,
                               const std::vector<std::complex<evalType> > &D1n,
                               const std::vector<std::complex<evalType> > &Zeta,
                               const std::vector<std::complex<evalType> > &D3n,
                               const std::vector<evalType> &Pi,
                               const std::vector<evalType> &Tau,
                               std::vector<std::complex<evalType> > &E,
                               std::vector<std::complex<evalType> > &H);

    bool isExpCoeffsCalc_ = false;
    bool isScaCoeffsCalc_ = false;
    bool isMieCalculated_ = false;

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
    /// Store result
    FloatType Qsca_ = 0.0, Qext_ = 0.0, Qabs_ = 0.0, Qbk_ = 0.0, Qpr_ = 0.0, asymmetry_factor_ = 0.0, albedo_ = 0.0;
    std::vector<std::vector< std::complex<FloatType> > > E_, H_;  // {X[], Y[], Z[]}
    std::vector<std::vector< std::complex<FloatType> > > Es_, Hs_;  // {X[], Y[], Z[]}
    std::vector<FloatType> Eabs_, Habs_;  // {X[], Y[], Z[]}
    std::vector<std::complex<FloatType> > S1_, S2_;
    void calcMieSeriesNeededToConverge(const FloatType Rho);
    void calcPiTauAllTheta(const double from_Theta,
                           const double to_Theta,
                           std::vector<std::vector<FloatType>> &Pi,
                           std::vector<std::vector<FloatType>> &Tau);
    void calcRadialOnlyDependantFunctions(const double from_Rho,
                                          const double to_Rho,
                                          const bool isIgnoreAvailableNmax,
                                          std::vector<std::vector<std::complex<FloatType>>> &Psi,
                                          std::vector<std::vector<std::complex<FloatType>>> &D1n,
                                          std::vector<std::vector<std::complex<FloatType>>> &Zeta,
                                          std::vector<std::vector<std::complex<FloatType>>> &D3n);
    void convertFieldsFromSphericalToCartesian();
  };  // end of class MultiLayerMie

}  // end of namespace nmie
#endif  // SRC_NMIE_HPP_
