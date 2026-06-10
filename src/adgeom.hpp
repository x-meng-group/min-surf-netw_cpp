// Small 3-vector helpers over autodiff scalars (AD), used when assembling the
// energy functional residuals.
#pragma once

#include "autodiff.hpp"
#include "linalg.hpp"

namespace msn {

struct Vec3AD {
    AD x, y, z;
    Vec3AD() = default;
    Vec3AD(AD x_, AD y_, AD z_) : x(x_), y(y_), z(z_) {}

    Vec3AD operator+(const Vec3AD& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3AD operator-(const Vec3AD& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3AD operator*(double s) const { return {x * s, y * s, z * s}; }
};

inline AD dot(const Vec3AD& a, const Vec3AD& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline AD normSq(const Vec3AD& a) { return adsquare(a.x) + adsquare(a.y) + adsquare(a.z); }

}  // namespace msn
