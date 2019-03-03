#include <pybind11/pybind11.h>
#include "pb11_nmie.hpp"

namespace py = pybind11;

// wrap as Python module
PYBIND11_MODULE(scattnlay_mp, m)
{
  m.doc() = "The Python version of scattnlay";

  m.def("scattcoeffs", &scattcoeffs,
        "Calculate the scattering coefficients, required to calculate both the near- and far-field parameters.",
        py::arg("x"), py::arg("m"), py::arg("nmax") = -1, py::arg("pl") = -1);

  m.def("scattnlay", &scattnlay,
        "Calculate the scattering parameters and amplitudes.",
        py::arg("x"), py::arg("m"), py::arg("theta") = py::array_t<double>(0), py::arg("nmax") = -1, py::arg("pl") = -1);

  m.def("fieldnlay", &fieldnlay,
        "Calculate the complex electric and magnetic field in the surroundings and inside the particle.",
        py::arg("x"), py::arg("m"), py::arg("coords"), py::arg("nmax") = -1, py::arg("pl") = -1);
}
