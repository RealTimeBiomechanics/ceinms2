#ifndef ceinms_CubicSpline_h
#define ceinms_CubicSpline_h
#include <vector>
#include "ceinms2/Curve.h"
#include "ceinms2/Types.h"

namespace ceinms {
/*This class wraps the more generic Curve class and makes it compatible with 
boost autodiff*/

class CubicSpline {
  public:
    struct Parameters {
        std::vector<DoubleT> x;
        std::vector<DoubleT> y;
    };
    CubicSpline() = default;
    CubicSpline(const Parameters &p) { set(p); }
    CubicSpline(const std::vector<DoubleT> &x, const std::vector<DoubleT> &y) { set({ x, y }); }

    constexpr DoubleT get(DoubleT x) const {
        x = std::min(x, spline_.getMaxX());
        x = std::max(x, spline_.getMinX());
        return spline_.getValue(x);
    }

    constexpr DoubleT getFirstDerivative(DoubleT x) const {
        x = std::min(x, spline_.getMaxX());
        x = std::max(x, spline_.getMinX());
        return spline_.getFirstDerivative(x);
    }

    constexpr DoubleT getSecondDerivative(DoubleT x) const {
        x = std::min(x, spline_.getMaxX());
        x = std::max(x, spline_.getMinX());
        return spline_.getSecondDerivative(x);
    }

    template<typename T>
    constexpr auto get(const T &cr) const {
        using root_type = typename T::root_type;
        constexpr size_t order = T::order_sum;
        root_type const d0 = get(static_cast<root_type>(cr));
        if constexpr (order == 0)
            return fvar(d0);
        else {
            const root_type d1 = getFirstDerivative(static_cast<root_type>(cr));
            const root_type d2 = getSecondDerivative(static_cast<root_type>(cr));
            const root_type derivatives[3]{ d0, d1, d2 };
            return cr.apply_derivatives(order, [&derivatives](size_t i) { return derivatives[i]; });
        }
    }

    void set(const Parameters &p) { spline_.resetPointsWith(p.x, p.y); }

  private:
    Curve<CurveMode::Offline> spline_;
};
}// namespace ceinms

#endif