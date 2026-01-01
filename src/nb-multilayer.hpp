#pragma once
#include "nb-helpers.hpp"
#include "nmie.hpp"

namespace nmie {
template <typename FloatType>
class PyMultiLayerMie : public MultiLayerMie<FloatType> {
public:
    // Verify: Call GetS1() from Python and check 's1.base' is not None (shows it shares memory)
    auto GetS1() { 
        return MoveVectorToNdarray(ConvertComplexVector<double>(this->S1_)); 
    }

    // Advanced: multi-dimensional field data
    auto GetFieldE() {
        return MoveVector2DToNdarray(ConvertComplexVectorVector<double>(this->E_));
    }

    void SetLayersSize(double x) {
        std::vector<FloatType> v = {static_cast<FloatType>(x)};
        this->MultiLayerMie<FloatType>::SetLayersSize(v);
    }

    void SetLayersSize(nb::ndarray<double, nb::c_contig> x) {
        auto v = NdarrayToVector(x);
        if constexpr (std::is_same_v<FloatType, double>) {
            this->MultiLayerMie<FloatType>::SetLayersSize(v);
        } else {
            std::vector<FloatType> v_conv;
            v_conv.reserve(v.size());
            for (const auto& val : v) {
                v_conv.push_back(static_cast<FloatType>(val));
            }
            this->MultiLayerMie<FloatType>::SetLayersSize(v_conv);
        }
    }

    void SetLayersIndex(std::complex<double> m) {
        std::vector<std::complex<FloatType>> v = {
            std::complex<FloatType>(static_cast<FloatType>(m.real()), static_cast<FloatType>(m.imag()))
        };
        this->MultiLayerMie<FloatType>::SetLayersIndex(v);
    }

    void SetLayersIndex(nb::ndarray<std::complex<double>, nb::c_contig> x) {
        auto v = NdarrayToVector(x);
        if constexpr (std::is_same_v<FloatType, double>) {
            this->MultiLayerMie<FloatType>::SetLayersIndex(v);
        } else {
            std::vector<std::complex<FloatType>> v_conv;
            v_conv.reserve(v.size());
            for (const auto& val : v) {
                v_conv.emplace_back(static_cast<FloatType>(val.real()), static_cast<FloatType>(val.imag()));
            }
            this->MultiLayerMie<FloatType>::SetLayersIndex(v_conv);
        }
    }

    void SetAngles(nb::ndarray<double, nb::c_contig> x) {
        auto v = NdarrayToVector(x);
        if constexpr (std::is_same_v<FloatType, double>) {
            this->MultiLayerMie<FloatType>::SetAngles(v);
        } else {
            std::vector<FloatType> v_conv;
            v_conv.reserve(v.size());
            for (const auto& val : v) {
                v_conv.push_back(static_cast<FloatType>(val));
            }
            this->MultiLayerMie<FloatType>::SetAngles(v_conv);
        }
    }

    void SetFieldCoords(nb::ndarray<double, nb::c_contig> xp, nb::ndarray<double, nb::c_contig> yp, nb::ndarray<double, nb::c_contig> zp) {
        auto v_xp = NdarrayToVector(xp);
        auto v_yp = NdarrayToVector(yp);
        auto v_zp = NdarrayToVector(zp);
        
        if constexpr (std::is_same_v<FloatType, double>) {
            this->MultiLayerMie<FloatType>::SetFieldCoords({v_xp, v_yp, v_zp});
        } else {
            std::vector<FloatType> v_xp_conv, v_yp_conv, v_zp_conv;
            v_xp_conv.reserve(v_xp.size());
            v_yp_conv.reserve(v_yp.size());
            v_zp_conv.reserve(v_zp.size());
            for (const auto& val : v_xp) v_xp_conv.push_back(static_cast<FloatType>(val));
            for (const auto& val : v_yp) v_yp_conv.push_back(static_cast<FloatType>(val));
            for (const auto& val : v_zp) v_zp_conv.push_back(static_cast<FloatType>(val));
            this->MultiLayerMie<FloatType>::SetFieldCoords({v_xp_conv, v_yp_conv, v_zp_conv});
        }
    }

    auto GetS2() { return MoveVectorToNdarray(ConvertComplexVector<double>(this->S2_)); }
    auto GetAn() { return MoveVectorToNdarray(ConvertComplexVector<double>(this->an_)); }
    auto GetBn() { return MoveVectorToNdarray(ConvertComplexVector<double>(this->bn_)); }
    
    auto GetFieldH() {
        return MoveVector2DToNdarray(ConvertComplexVectorVector<double>(this->H_));
    }

    auto GetFieldEabs() { return MoveVectorToNdarray(ConvertVector<double>(this->Eabs_)); }
    auto GetFieldHabs() { return MoveVectorToNdarray(ConvertVector<double>(this->Habs_)); }

    auto GetLayerAn() { return MoveVector2DToNdarray(ConvertComplexVectorVector<double>(this->aln_)); }
    auto GetLayerBn() { return MoveVector2DToNdarray(ConvertComplexVectorVector<double>(this->bln_)); }
    auto GetLayerCn() { return MoveVector2DToNdarray(ConvertComplexVectorVector<double>(this->cln_)); }
    auto GetLayerDn() { return MoveVector2DToNdarray(ConvertComplexVectorVector<double>(this->dln_)); }

    // Scalar getters (simple wrappers)
    double GetQext() { return static_cast<double>(this->Qext_); }
    double GetQsca() { return static_cast<double>(this->Qsca_); }
    double GetQabs() { return static_cast<double>(this->Qabs_); }
    double GetQbk() { return static_cast<double>(this->Qbk_); }
    double GetQpr() { return static_cast<double>(this->Qpr_); }
    double GetAsymmetryFactor() { return static_cast<double>(this->asymmetry_factor_); }
    double GetAlbedo() { return static_cast<double>(this->albedo_); }
};
}
