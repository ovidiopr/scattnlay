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

#include <string>
#include "nmie-basic.hpp"
#include "nmie-nearfield.hpp"
#include "nmie.hpp"
#include "pb11-multilayer.hpp"

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

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
