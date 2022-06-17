#include <string>
#include "nmie-basic.hpp"
#include "nmie-nearfield.hpp"
#include "nmie.hpp"
#include "pb11-helpers.hpp"

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace nmie {
template <typename FloatType>
class PyMultiLayerMie : public MultiLayerMie<FloatType> {
 public:
  template <typename outputType>
  py::array_t<std::complex<outputType>> GetS1();
  template <typename outputType>
  py::array_t<std::complex<outputType>> GetS2();
  template <typename outputType>
  py::array_t<std::complex<outputType>> GetAn();
  template <typename outputType>
  py::array_t<std::complex<outputType>> GetBn();
  template <typename outputType>
  py::array GetLayerAn();
  template <typename outputType>
  py::array GetLayerBn();
  template <typename outputType>
  py::array GetLayerCn();
  template <typename outputType>
  py::array GetLayerDn();
  template <typename outputType>
  py::array GetFieldE();
  template <typename outputType>
  py::array GetFieldH();
  template <typename outputType>
  py::array_t<outputType> GetFieldEabs();
  template <typename outputType>
  py::array_t<outputType> GetFieldHabs();

  template <typename inputType>
  void SetLayersSize(
      const py::array_t<inputType, py::array::c_style | py::array::forcecast>&
          py_layer_size);
  template <typename inputType>
  void SetLayersIndex(
      const py::array_t<std::complex<inputType>,
                        py::array::c_style | py::array::forcecast>& py_index);
  template <typename inputType>
  void SetAngles(
      const py::array_t<inputType, py::array::c_style | py::array::forcecast>&
          py_angles);
  void SetFieldCoords(
      const py::array_t<double, py::array::c_style | py::array::forcecast>&
          py_Xp,
      const py::array_t<double, py::array::c_style | py::array::forcecast>&
          py_Yp,
      const py::array_t<double, py::array::c_style | py::array::forcecast>&
          py_Zp);
};

// Python interface
template <typename FloatType>
template <typename inputType>
void PyMultiLayerMie<FloatType>::SetLayersSize(
    const py::array_t<inputType, py::array::c_style | py::array::forcecast>&
        py_layer_size) {
  auto layer_size_dp = Py2Vector<inputType>(py_layer_size);
  this->MultiLayerMie<FloatType>::SetLayersSize(
      ConvertVector<FloatType>(layer_size_dp));
}

template <typename FloatType>
template <typename inputType>
void PyMultiLayerMie<FloatType>::SetLayersIndex(
    const py::array_t<std::complex<inputType>,
                      py::array::c_style | py::array::forcecast>& py_index) {
  auto index_dp = Py2Vector<std::complex<inputType>>(py_index);
  this->MultiLayerMie<FloatType>::SetLayersIndex(
      ConvertComplexVector<FloatType>(index_dp));
}

template <typename FloatType>
template <typename inputType>
void PyMultiLayerMie<FloatType>::SetAngles(
    const py::array_t<inputType, py::array::c_style | py::array::forcecast>&
        py_angles) {
  auto angles_dp = Py2Vector<inputType>(py_angles);
  this->MultiLayerMie<FloatType>::SetAngles(
      ConvertVector<FloatType>(angles_dp));
}

template <typename FloatType>
template <typename outputType>
py::array_t<std::complex<outputType>> PyMultiLayerMie<FloatType>::GetS1() {
  return VectorComplex2Py<FloatType, outputType>(
      this->MultiLayerMie<FloatType>::GetS1());
}

template <typename FloatType>
template <typename outputType>
py::array_t<std::complex<outputType>> PyMultiLayerMie<FloatType>::GetS2() {
  return VectorComplex2Py<FloatType, outputType>(
      this->MultiLayerMie<FloatType>::GetS2());
}

template <typename FloatType>
template <typename outputType>
py::array_t<outputType> PyMultiLayerMie<FloatType>::GetFieldEabs() {
  return Vector2Py(
      ConvertVector<double>(this->MultiLayerMie<FloatType>::GetFieldEabs()));
}

template <typename FloatType>
template <typename outputType>
py::array_t<outputType> PyMultiLayerMie<FloatType>::GetFieldHabs() {
  return Vector2Py(
      ConvertVector<double>(this->MultiLayerMie<FloatType>::GetFieldHabs()));
}

template <typename FloatType>
template <typename outputType>
py::array_t<std::complex<outputType>> PyMultiLayerMie<FloatType>::GetAn() {
  return VectorComplex2Py<FloatType, outputType>(
      this->MultiLayerMie<FloatType>::GetAn());
}

template <typename FloatType>
template <typename outputType>
py::array_t<std::complex<outputType>> PyMultiLayerMie<FloatType>::GetBn() {
  return VectorComplex2Py<FloatType, outputType>(
      this->MultiLayerMie<FloatType>::GetBn());
}

template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetFieldE() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetFieldE()));
}

template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetFieldH() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetFieldH()));
}

template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetLayerAn() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetLayerAn()));
}

template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetLayerBn() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetLayerBn()));
}

template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetLayerCn() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetLayerCn()));
}

template <typename FloatType>
template <typename outputType>
py::array PyMultiLayerMie<FloatType>::GetLayerDn() {
  return Vector2DComplex2Py<std::complex<outputType>>(
      ConvertComplexVectorVector<outputType>(
          this->MultiLayerMie<FloatType>::GetLayerDn()));
}

template <typename FloatType>
void PyMultiLayerMie<FloatType>::SetFieldCoords(
    const py::array_t<double, py::array::c_style | py::array::forcecast>& py_Xp,
    const py::array_t<double, py::array::c_style | py::array::forcecast>& py_Yp,
    const py::array_t<double, py::array::c_style | py::array::forcecast>&
        py_Zp) {
  auto c_Xp = Py2Vector<double>(py_Xp);
  auto c_Yp = Py2Vector<double>(py_Yp);
  auto c_Zp = Py2Vector<double>(py_Zp);
  this->MultiLayerMie<FloatType>::SetFieldCoords(
      {ConvertVector<FloatType>(c_Xp), ConvertVector<FloatType>(c_Yp),
       ConvertVector<FloatType>(c_Zp)});
}
}  // namespace nmie

template <typename T>
void declare_nmie(py::module& m, const std::string& typestr) {
  using mie_typed = nmie::PyMultiLayerMie<nmie::FloatType>;
  std::string pyclass_name = std::string("mie") + typestr;
  py::class_<mie_typed>(m, pyclass_name.c_str(), py::buffer_protocol(),
                        py::dynamic_attr())
      .def(py::init<>())
      .def("GetPECLayer", &mie_typed::GetPECLayer)
      .def("SetPECLayer", &mie_typed::SetPECLayer)
      .def("SetMaxTerms", &mie_typed::SetMaxTerms)
      .def("GetMaxTerms", &mie_typed::GetMaxTerms)
      .def("SetModeNmaxAndType", &mie_typed::SetModeNmaxAndType)
      .def("RunMieCalculation", &mie_typed::RunMieCalculation)
      .def("calcScattCoeffs", &mie_typed::calcScattCoeffs)
      .def("calcExpanCoeffs", &mie_typed::calcExpanCoeffs)
      .def("RunFieldCalculation", &mie_typed::RunFieldCalculation,
           py::arg("isMarkUnconverged") = true)
      .def("RunFieldCalculationPolar", &mie_typed::RunFieldCalculationPolar,
           py::arg("outer_arc_points") = 1, py::arg("radius_points") = 1,
           py::arg("from_Rho") = 0, py::arg("to_Rho") = 1,
           py::arg("from_Theta") = 0, py::arg("to_Theta") = 3.14159265358979323,
           py::arg("from_Phi") = 0, py::arg("to_Phi") = 3.14159265358979323,
           py::arg("isMarkUnconverged") = true, py::arg("nmax_in") = -1)
      .def("GetPECLayer", &mie_typed::GetPECLayer)

      .def("SetLayersSize",
           static_cast<void (mie_typed::*)(
               const py::array_t<double,
                                 py::array::c_style | py::array::forcecast>&)>(
               &mie_typed::SetLayersSize))
      .def("SetLayersIndex",
           static_cast<void (mie_typed::*)(
               const py::array_t<std::complex<double>,
                                 py::array::c_style | py::array::forcecast>&)>(
               &mie_typed::SetLayersIndex))
      .def("SetAngles",
           static_cast<void (mie_typed::*)(
               const py::array_t<double,
                                 py::array::c_style | py::array::forcecast>&)>(
               &mie_typed::SetAngles))
      .def("SetFieldCoords",
           static_cast<void (mie_typed::*)(
               const py::array_t<double, py::array::c_style |
                                             py::array::forcecast>& py_Xp,
               const py::array_t<double, py::array::c_style |
                                             py::array::forcecast>& py_Yp,
               const py::array_t<double, py::array::c_style |
                                             py::array::forcecast>& py_Zp)>(
               &mie_typed::SetFieldCoords))
      .def("GetQext", &mie_typed::GetQext<double>)
      .def("GetQsca", &mie_typed::GetQsca<double>)
      .def("GetQabs", &mie_typed::GetQabs<double>)
      .def("GetQpr", &mie_typed::GetQpr<double>)
      .def("GetQbk", &mie_typed::GetQbk<double>)
      .def("GetAsymmetryFactor", &mie_typed::GetAsymmetryFactor<double>)
      .def("GetAlbedo", &mie_typed::GetAlbedo<double>)
      .def("GetS1", &mie_typed::GetS1<double>)
      .def("GetS2", &mie_typed::GetS2<double>)
      .def("GetAn", &mie_typed::GetAn<double>)
      .def("GetBn", &mie_typed::GetBn<double>)
      .def("GetFieldE", &mie_typed::GetFieldE<double>)
      .def("GetFieldH", &mie_typed::GetFieldH<double>)
      .def("GetFieldEabs", &mie_typed::GetFieldEabs<double>)
      .def("GetFieldHabs", &mie_typed::GetFieldHabs<double>)
      .def("GetFieldConvergence", &mie_typed::GetFieldConvergence)
      .def("GetLayerAn", &mie_typed::GetLayerAn<double>)
      .def("GetLayerBn", &mie_typed::GetLayerBn<double>)
      .def("GetLayerCn", &mie_typed::GetLayerCn<double>)
      .def("GetLayerDn", &mie_typed::GetLayerDn<double>);
}

// wrap as Python module
#ifdef MULTI_PRECISION
std::string precision_name = "_mp";
PYBIND11_MODULE(scattnlay_mp, m)
#else
std::string precision_name = "_dp";
PYBIND11_MODULE(scattnlay_dp, m)
#endif  // MULTI_PRECISION
{
  m.doc() = "The Python version of scattnlay";
  declare_nmie<nmie::FloatType>(m, precision_name);
}
