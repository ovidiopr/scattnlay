#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>

#include "nmie.hpp"
#include "nmie-basic.hpp"
#include "nmie-nearfield.hpp"
#include "mesomie.hpp"
#include "nb-multilayer.hpp"

namespace nb = nanobind;

template <typename FloatType, typename Engine>
void declare_mie(nb::module_& m, const std::string& pyclass_name) {
    using mie_typed = nmie::PyMultiLayerMie<FloatType, Engine>;
    
    nb::class_<mie_typed>(m, pyclass_name.c_str())
        .def(nb::init<>())
        .def("SetLayersSize", nb::overload_cast<double>(&mie_typed::SetLayersSize))
        .def("SetLayersSize", nb::overload_cast<nb::ndarray<double, nb::c_contig>>(&mie_typed::SetLayersSize))
        .def("SetLayersIndex", nb::overload_cast<std::complex<double>>(&mie_typed::SetLayersIndex))
        .def("SetLayersIndex", nb::overload_cast<nb::ndarray<std::complex<double>, nb::c_contig>>(&mie_typed::SetLayersIndex))
        .def("SetAngles", &mie_typed::SetAngles)
        .def("SetFieldCoords", &mie_typed::SetFieldCoords)
        .def("SetPECLayer", &mie_typed::SetPECLayer)
        .def("GetPECLayer", &mie_typed::GetPECLayer)
        .def("SetMaxTerms", &mie_typed::SetMaxTerms)
        .def("GetMaxTerms", &mie_typed::GetMaxTerms)
        .def("SetModeNmaxAndType", &mie_typed::SetModeNmaxAndType)
        .def("RunMieCalculation", &mie_typed::RunMieCalculation)
        .def("RunFieldCalculation", &mie_typed::RunFieldCalculation, nb::arg("isMarkUnconverged") = true)
        .def("RunFieldCalculationPolar", &mie_typed::RunFieldCalculationPolar,
             nb::arg("outer_arc_points") = 1, nb::arg("radius_points") = 1,
             nb::arg("from_Rho") = 0, nb::arg("to_Rho") = 1,
             nb::arg("from_Theta") = 0, nb::arg("to_Theta") = 3.14159265358979323,
             nb::arg("from_Phi") = 0, nb::arg("to_Phi") = 3.14159265358979323,
             nb::arg("isMarkUnconverged") = true, nb::arg("nmax_in") = -1)
        .def("calcScattCoeffs", &mie_typed::calcScattCoeffs)
        .def("calcExpanCoeffs", &mie_typed::calcExpanCoeffs)
        .def("GetQext", &mie_typed::GetQext)
        .def("GetQsca", &mie_typed::GetQsca)
        .def("GetQabs", &mie_typed::GetQabs)
        .def("GetQbk", &mie_typed::GetQbk)
        .def("GetQpr", &mie_typed::GetQpr)
        .def("GetAsymmetryFactor", &mie_typed::GetAsymmetryFactor)
        .def("GetAlbedo", &mie_typed::GetAlbedo)
        .def("GetS1", &mie_typed::GetS1)
        .def("GetS2", &mie_typed::GetS2)
        .def("GetAn", &mie_typed::GetAn)
        .def("GetBn", &mie_typed::GetBn)
        .def("GetFieldE", &mie_typed::GetFieldE)
        .def("GetFieldH", &mie_typed::GetFieldH)
        .def("GetFieldEabs", &mie_typed::GetFieldEabs)
        .def("GetFieldHabs", &mie_typed::GetFieldHabs)
        .def("GetFieldConvergence", &mie_typed::GetFieldConvergence)
        .def("GetLayerAn", &mie_typed::GetLayerAn)
        .def("GetLayerBn", &mie_typed::GetLayerBn)
        .def("GetLayerCn", &mie_typed::GetLayerCn)
        .def("GetLayerDn", &mie_typed::GetLayerDn);
}

template <typename FloatType>
void declare_mesomie(nb::module_& m, const std::string& pyclass_name) {
    using mesomie = nmie::MesoMie<FloatType>;
    nb::class_<mesomie>(m, pyclass_name.c_str())
        .def(nb::init<>())
        .def("calc_Q", &mesomie::calc_Q)
        .def("calc_ab", &mesomie::calc_ab, nb::arg("R") = 1, nb::arg("xd") = 1,
             nb::arg("xm") = 1, nb::arg("eps_d") = 1, nb::arg("eps_m") = 1,
             nb::arg("d_parallel") = 0, nb::arg("d_perp") = 0)
        .def("GetQext", &mesomie::template GetQext<double>)
        .def("GetQsca", &mesomie::template GetQsca<double>);
}

#ifdef MULTI_PRECISION
std::string precision_name = "_mp";
NB_MODULE(scattnlay_mp, m) {
    m.doc() = "The Python version of scattnlay (nanobind)";
    declare_mie<nmie::FloatType, nmie::DefaultEngine<nmie::FloatType>>(m, "mie_mp");
    declare_mesomie<nmie::FloatType>(m, "mesomie_mp");
}
#else
std::string precision_name = "_dp";
NB_MODULE(scattnlay_dp, m) {
    m.doc() = "The Python version of scattnlay (nanobind)";
    
    declare_mie<double, nmie::DefaultEngine<double>>(m, "mie_dp");
    declare_mesomie<double>(m, "mesomie_dp");

#ifdef WITH_HWY
    declare_mie<double, nmie::ScalarEngine<double>>(m, "mie_scalar");
#endif
}
#endif

#ifdef BUILD_SIMD_MODULE
#include "nmie-batch.hpp"

// High-performance batch processor
nb::dict RunMieBatchPy(
    nb::ndarray<double, nb::numpy, nb::c_contig, nb::device::cpu> x, 
    nb::ndarray<std::complex<double>, nb::numpy, nb::c_contig, nb::device::cpu> m,
    nb::ndarray<double, nb::numpy, nb::c_contig, nb::device::cpu> theta
) {
    if (x.size() != m.size()) {
        throw std::invalid_argument("Size parameters (x) and indices (m) must have the same length.");
    }

    nmie::MieBatchInput input;
    // Fast assignment from NumPy pointers to C++ vectors
    input.x.assign(x.data(), x.data() + x.size());
    input.m.assign(m.data(), m.data() + m.size());

    if (theta.size() > 0) {
        input.theta.assign(theta.data(), theta.data() + theta.size());
    }

    // Heavy lifting: Google Highway SIMD execution
    auto output = nmie::RunMieBatch<double>(input);

    // Build return dictionary using Zero-Copy "Move" helpers
    nb::dict res;
    res["Qext"]   = MoveVectorToNdarray(std::move(output.Qext));
    res["Qsca"]   = MoveVectorToNdarray(std::move(output.Qsca));
    res["Qabs"]   = MoveVectorToNdarray(std::move(output.Qabs));
    res["Qbk"]    = MoveVectorToNdarray(std::move(output.Qbk));
    res["Qpr"]    = MoveVectorToNdarray(std::move(output.Qpr));
    res["g"]      = MoveVectorToNdarray(std::move(output.g));
    res["Albedo"] = MoveVectorToNdarray(std::move(output.Albedo));
    
    if (!output.S1.empty()) {
        res["S1"] = MoveVector2DToNdarray(std::move(output.S1));
        res["S2"] = MoveVector2DToNdarray(std::move(output.S2));
    }

    return res;
}

NB_MODULE(scattnlay_simd, m) {
    m.doc() = "Google Highway SIMD accelerated Mie solver (nanobind)";
    m.def("RunMieBatch", &RunMieBatchPy, 
          nb::arg("x"), nb::arg("m"), nb::arg("theta") = nb::none(),
          "Run Mie calculation for a batch of particles using SIMD");
}
#endif
