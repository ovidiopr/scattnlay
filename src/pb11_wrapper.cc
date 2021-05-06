#include <string>
#include <pybind11/pybind11.h>
#include "nmie-pybind11.hpp"
#include "nmie.hpp"

namespace py = pybind11;

template<typename T>
void declare_nmie(py::module &m, std::string &typestr) {
  using mie_typed = nmie::MultiLayerMie<nmie::FloatType>;
  std::string pyclass_name = std::string("mie") + typestr;
  py::class_<mie_typed>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
      .def(py::init<>())
      .def("GetPECLayer", &mie_typed::GetPECLayer)
//      .def("width",     &mie_typed::width)
//      .def("height",    &mie_typed::height)
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

  m.def("scattcoeffs", &py_ScattCoeffs,
        "Calculate the scattering coefficients, required to calculate both the near- and far-field parameters.",
        py::arg("x"), py::arg("m"), py::arg("nmax") = -1, py::arg("pl") = -1);

  m.def("expancoeffs", &py_ExpanCoeffs,
        "Calculate the expansion coefficients, required to calculate the near-field parameters.",
        py::arg("x"), py::arg("m"), py::arg("nmax") = -1, py::arg("pl") = -1);

  m.def("scattnlay", &py_scattnlay,
        "Calculate the scattering parameters and amplitudes.",
        py::arg("x"), py::arg("m"), py::arg("theta") = py::array_t<double>(0), py::arg("nmax") = -1, py::arg("pl") = -1);

  m.def("fieldnlay", &py_fieldnlay,
        "Calculate the complex electric and magnetic field in the surroundings and inside the particle.",
        py::arg("x"), py::arg("m"), py::arg("xp"), py::arg("yp"), py::arg("zp"), py::arg("nmax") = -1, py::arg("pl") = -1);
}
