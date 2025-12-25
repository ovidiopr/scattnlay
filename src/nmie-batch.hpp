#ifndef SRC_NMIE_BATCH_HPP_
#define SRC_NMIE_BATCH_HPP_

#include <vector>
#include <complex>
#include <algorithm>
#include "nmie-precision.hpp"
#include "special-functions-impl.hpp"

namespace nmie {

struct MieBatchInput {
  std::vector<double> x;
  std::vector<std::complex<double>> m;
};

struct MieBatchOutput {
  std::vector<double> Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
};

// Dispatcher function template
template <typename FloatType>
MieBatchOutput RunMieBatch(const MieBatchInput& input);

} // namespace nmie

#endif // SRC_NMIE_BATCH_HPP_

// Highway dynamic dispatch boilerplate
#ifndef HWY_TARGET_TOGGLE
#include <hwy/highway.h>
#include <hwy/contrib/math/math-inl.h>
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/nmie-batch.hpp"
#include <hwy/foreach_target.h>  // IWYU pragma: keep

#ifdef HWY_TARGET_TOGGLE
#include "nmie-batch-inl.h"
#endif

#if HWY_ONCE
namespace nmie {
using namespace hwy;
HWY_EXPORT(RunMieBatchDouble);

template <>
inline MieBatchOutput RunMieBatch<double>(const MieBatchInput& input) {
    return HWY_DYNAMIC_DISPATCH(RunMieBatchDouble)(input);
}

} // namespace nmie
#endif
