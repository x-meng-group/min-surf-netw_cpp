// Order-1 (piecewise-linear) interpolation of a 3D path, replacing
// Mathematica's Interpolation[..., InterpolationOrder -> 1]. Supports
// evaluation outside the knot range via linear extrapolation from the end
// segments (Mathematica's default behaviour for order-1 interpolation).
#pragma once

#include <vector>

#include "linalg.hpp"

namespace msn {

class LinearInterp {
public:
    LinearInterp() = default;
    LinearInterp(std::vector<double> knots, std::vector<Vec3> values)
        : t_(std::move(knots)), v_(std::move(values)) {}

    // Evaluate at parameter t (extrapolates beyond endpoints).
    Vec3 operator()(double t) const {
        const std::size_t n = t_.size();
        if (n == 0) return {0, 0, 0};
        if (n == 1) return v_[0];
        std::size_t i = 0;
        if (t <= t_[0]) {
            i = 0;
        } else if (t >= t_[n - 1]) {
            i = n - 2;
        } else {
            // find segment [t_[i], t_[i+1]] containing t
            while (i + 1 < n && t_[i + 1] < t) ++i;
            if (i > n - 2) i = n - 2;
        }
        double dt = t_[i + 1] - t_[i];
        double a = (dt != 0.0) ? (t - t_[i]) / dt : 0.0;
        return v_[i] + (v_[i + 1] - v_[i]) * a;
    }

    double tMin() const { return t_.front(); }
    double tMax() const { return t_.back(); }

private:
    std::vector<double> t_;
    std::vector<Vec3> v_;
};

// Linearly resample the polyline (p0..p_{m-1}) so the result has `nOut` evenly
// spaced samples in arclength-parameter space, matching ArrayResample of a
// 2-point segment used in the original code's boundary generation.
inline std::vector<Vec3> resampleSegment(const Vec3& p0, const Vec3& p1, int nOut) {
    std::vector<Vec3> out;
    out.reserve(nOut);
    if (nOut <= 1) {
        out.push_back(p0);
        return out;
    }
    for (int i = 0; i < nOut; ++i) {
        double a = static_cast<double>(i) / (nOut - 1);
        out.push_back(p0 + (p1 - p0) * a);
    }
    return out;
}

}  // namespace msn
