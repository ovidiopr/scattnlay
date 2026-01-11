#pragma once
#include "nb-helpers.hpp"
#include "nmie.hpp"

namespace nmie {
template <typename FloatType, typename Engine = DefaultEngine<FloatType>>
class PyMultiLayerMie : public MultiLayerMie<FloatType, Engine> {
public:
    // Verify: Call GetS1() from Python and check 's1.base' is not None (shows it shares memory)
    auto GetS1() { 
        if constexpr (std::is_same_v<FloatType, double>) {
            size_t shape[1] = { this->S1_.size() };
            return nb::ndarray<nb::numpy, std::complex<double>, nb::shape<-1>>(
                this->S1_.data(), 1, shape, nb::cast(this)
            );
        } else {
            return MoveVectorToNdarray(ConvertComplexVector<double>(this->S1_)); 
        }
    }

    // Advanced: multi-dimensional field data
    auto GetFieldE() {
        return MoveVector2DToNdarray(ConvertComplexVectorVector<double>(this->E_));
    }

    void SetLayersSize(double x) {
        std::vector<FloatType> v = {static_cast<FloatType>(x)};
        this->MultiLayerMie<FloatType, Engine>::SetLayersSize(v);
    }

    void SetLayersSize(nb::ndarray<double, nb::c_contig> x) {
        auto v = NdarrayToVector(x);
        if constexpr (std::is_same_v<FloatType, double>) {
            this->MultiLayerMie<FloatType, Engine>::SetLayersSize(v);
        } else {
            std::vector<FloatType> v_conv;
            v_conv.reserve(v.size());
            for (const auto& val : v) {
                v_conv.push_back(static_cast<FloatType>(val));
            }
            this->MultiLayerMie<FloatType, Engine>::SetLayersSize(v_conv);
        }
    }

    void SetLayersIndex(std::complex<double> m) {
        std::vector<std::complex<FloatType>> v = {
            std::complex<FloatType>(static_cast<FloatType>(m.real()), static_cast<FloatType>(m.imag()))
        };
        this->MultiLayerMie<FloatType, Engine>::SetLayersIndex(v);
    }

    void SetLayersIndex(nb::ndarray<std::complex<double>, nb::c_contig> x) {
        auto v = NdarrayToVector(x);
        if constexpr (std::is_same_v<FloatType, double>) {
            this->MultiLayerMie<FloatType, Engine>::SetLayersIndex(v);
        } else {
            std::vector<std::complex<FloatType>> v_conv;
            v_conv.reserve(v.size());
            for (const auto& val : v) {
                v_conv.emplace_back(static_cast<FloatType>(val.real()), static_cast<FloatType>(val.imag()));
            }
            this->MultiLayerMie<FloatType, Engine>::SetLayersIndex(v_conv);
        }
    }

    void SetAngles(nb::ndarray<double, nb::c_contig> x) {
        auto v = NdarrayToVector(x);
        if constexpr (std::is_same_v<FloatType, double>) {
            this->MultiLayerMie<FloatType, Engine>::SetAngles(v);
        } else {
            std::vector<FloatType> v_conv;
            v_conv.reserve(v.size());
            for (const auto& val : v) {
                v_conv.push_back(static_cast<FloatType>(val));
            }
            this->MultiLayerMie<FloatType, Engine>::SetAngles(v_conv);
        }
    }

    void SetFieldCoords(nb::ndarray<double, nb::c_contig> xp, nb::ndarray<double, nb::c_contig> yp, nb::ndarray<double, nb::c_contig> zp) {
        // Constructor-based copy from pointer is fast, but we only want to do it once
        std::vector<FloatType> v_xp(xp.data(), xp.data() + xp.size());
        std::vector<FloatType> v_yp(yp.data(), yp.data() + yp.size());
        std::vector<FloatType> v_zp(zp.data(), zp.data() + zp.size());
        
        if constexpr (std::is_same_v<FloatType, double>) {
            this->MultiLayerMie<FloatType, Engine>::SetFieldCoords({std::move(v_xp), std::move(v_yp), std::move(v_zp)});
        } else {
            // For non-double types, we still need to convert, so we can't avoid the copy/conversion
            // But we can still use the new signature of SetFieldCoords
            this->MultiLayerMie<FloatType, Engine>::SetFieldCoords({std::move(v_xp), std::move(v_yp), std::move(v_zp)});
        }
    }

    auto GetS2() { 
        if constexpr (std::is_same_v<FloatType, double>) {
            size_t shape[1] = { this->S2_.size() };
            return nb::ndarray<nb::numpy, std::complex<double>, nb::shape<-1>>(
                this->S2_.data(), 1, shape, nb::cast(this)
            );
        } else {
            return MoveVectorToNdarray(ConvertComplexVector<double>(this->S2_)); 
        }
    }
    auto GetAn() { return MoveVectorToNdarray(ConvertComplexVector<double>(this->an_)); }
    auto GetBn() { return MoveVectorToNdarray(ConvertComplexVector<double>(this->bn_)); }
    
    auto GetFieldH() {
        return MoveVector2DToNdarray(ConvertComplexVectorVector<double>(this->H_));
    }

    auto GetFieldEabs() { 
        if constexpr (std::is_same_v<FloatType, double>) {
            size_t shape[1] = { this->Eabs_.size() };
            return nb::ndarray<nb::numpy, double, nb::shape<-1>>(
                this->Eabs_.data(), 1, shape, nb::cast(this)
            );
        } else {
            return MoveVectorToNdarray(ConvertVector<double>(this->Eabs_)); 
        }
    }
    auto GetFieldHabs() { 
        if constexpr (std::is_same_v<FloatType, double>) {
            size_t shape[1] = { this->Habs_.size() };
            return nb::ndarray<nb::numpy, double, nb::shape<-1>>(
                this->Habs_.data(), 1, shape, nb::cast(this)
            );
        } else {
            return MoveVectorToNdarray(ConvertVector<double>(this->Habs_)); 
        }
    }

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