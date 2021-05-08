#include <string>
#include "nmie-pybind11.hpp"
#include "nmie.hpp"
#include "nmie-basic.hpp"
#include "nmie-nearfield.hpp"

//py::class_<Pet>(m, "Pet")
//.def(py::init<const std::string &, int>())
//.def("set", static_cast<void (Pet::*)(int)>(&Pet::set), "Set the pet's age")
//.def("set", static_cast<void (Pet::*)(const std::string &)>(&Pet::set), "Set the pet's name");
//

template<typename T>
void declare_nmie(py::module &m, std::string &typestr) {
  using mie_typed = nmie::MultiLayerMie<nmie::FloatType>;
  std::string pyclass_name = std::string("mie") + typestr;
  py::class_<mie_typed>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
      .def(py::init<>())
      .def("GetPECLayer", &mie_typed::GetPECLayer)
      .def("SetPECLayer", &mie_typed::SetPECLayer)
      .def("SetMaxTerms", &mie_typed::SetMaxTerms)
      .def("GetMaxTerms", &mie_typed::GetMaxTerms)
      .def("SetModeNmaxAndType", &mie_typed::SetModeNmaxAndType)
      .def("RunMieCalculation", &mie_typed::RunMieCalculation)
      .def("calcScattCoeffs", &mie_typed::calcScattCoeffs)
      .def("calcExpanCoeffs", &mie_typed::calcExpanCoeffs)
      .def("RunFieldCalculation", &mie_typed::RunFieldCalculation)
      .def("RunFieldCalculationPolar", &mie_typed::RunFieldCalculationPolar)
      .def("GetPECLayer", &mie_typed::GetPECLayer)

      .def("SetLayersSize", static_cast<
               void (mie_typed::*)
                   (const py::array_t<double, py::array::c_style | py::array::forcecast>&)>
           (&mie_typed::SetLayersSize)
      )
      .def("SetLayersIndex", static_cast<
               void (mie_typed::*)
                   (const py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> &)>
           (&mie_typed::SetLayersIndex)
      )
      .def("SetAngles", static_cast<
               void (mie_typed::*)
                   (const py::array_t<double, py::array::c_style | py::array::forcecast>&)>
           (&mie_typed::SetAngles)
      )
      .def("SetFieldCoords", static_cast<
               void (mie_typed::*)
                   (
                       const py::array_t<double, py::array::c_style | py::array::forcecast> &py_Xp,
                       const py::array_t<double, py::array::c_style | py::array::forcecast> &py_Yp,
                       const py::array_t<double, py::array::c_style | py::array::forcecast> &py_Zp
                   )>
           (&mie_typed::SetFieldCoords)
      )
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
      ;
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

  m.def("expancoeffs", &py_ExpanCoeffs,
        "Calculate the expansion coefficients, required to calculate the near-field parameters.",
        py::arg("x"), py::arg("m"), py::arg("nmax") = -1, py::arg("pl") = -1);
}

