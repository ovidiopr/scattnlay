#include <pybind11/pybind11.h>
#include "nmie-pybind11.hpp"

namespace py = pybind11;

// wrap as Python module
#ifdef MULTI_PRECISION
PYBIND11_MODULE(scattnlay_mp_, m)
#else
PYBIND11_MODULE(scattnlay_, m)
#endif  // MULTI_PRECISION
{
  m.doc() = "The Python version of scattnlay";

  m.def("scattcoeffs_", &py_ScattCoeffs,
        "Calculate the scattering coefficients, required to calculate both the near- and far-field parameters.",
        py::arg("x"), py::arg("m"), py::arg("nmax") = -1, py::arg("pl") = -1);

  m.def("scattnlay_", &py_scattnlay,
        "Calculate the scattering parameters and amplitudes.",
        py::arg("x"), py::arg("m"), py::arg("theta") = py::array_t<double>(0), py::arg("nmax") = -1, py::arg("pl") = -1);

  m.def("fieldnlay_", &py_fieldnlay,
        "Calculate the complex electric and magnetic field in the surroundings and inside the particle.",
        py::arg("x"), py::arg("m"), py::arg("xp"), py::arg("yp"), py::arg("zp"), py::arg("nmax") = -1, py::arg("pl") = -1);
}
