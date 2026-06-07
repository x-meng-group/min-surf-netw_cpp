// Small fixed-size linear algebra utilities (Vec3, Mat3) used throughout the
// network solver. Header-only, no external dependencies.
#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

namespace msn {

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
    double& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator-() const { return {-x, -y, -z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
};

inline Vec3 operator*(double s, const Vec3& v) { return {v.x * s, v.y * s, v.z * s}; }

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalize(const Vec3& a) {
    double n = norm(a);
    if (n == 0.0) return {0.0, 0.0, 0.0};
    return a / n;
}

// 3x3 matrix stored row-major.
struct Mat3 {
    std::array<std::array<double, 3>, 3> m{};

    static Mat3 identity() {
        Mat3 r;
        r.m[0][0] = r.m[1][1] = r.m[2][2] = 1.0;
        return r;
    }

    Vec3 operator*(const Vec3& v) const {
        return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
    }

    Mat3 operator*(const Mat3& o) const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k) s += m[i][k] * o.m[k][j];
                r.m[i][j] = s;
            }
        return r;
    }

    Mat3 transpose() const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[j][i];
        return r;
    }
};

// Angle between two vectors (radians), matching Mathematica's VectorAngle.
inline double vectorAngle(const Vec3& a, const Vec3& b) {
    double na = norm(a), nb = norm(b);
    if (na == 0.0 || nb == 0.0) return 0.0;
    double c = dot(a, b) / (na * nb);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    return std::acos(c);
}

// Rotation matrix taking unit-ish vector `from` to `to` (rotation in the plane
// spanned by the two vectors). Mirrors Mathematica's RotationMatrix[{from,to}].
inline Mat3 rotationMatrixBetween(const Vec3& from, const Vec3& to) {
    Vec3 a = normalize(from);
    Vec3 b = normalize(to);
    Vec3 v = cross(a, b);
    double s = norm(v);
    double c = dot(a, b);
    if (s < 1e-15) {
        if (c > 0.0) return Mat3::identity();
        // 180-degree rotation: pick an axis orthogonal to a.
        Vec3 axis = std::abs(a.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
        axis = normalize(cross(a, axis));
        // Rodrigues with theta = pi: R = 2 axis axis^T - I.
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                r.m[i][j] = 2.0 * axis[i] * axis[j] - (i == j ? 1.0 : 0.0);
        return r;
    }
    Vec3 k = v / s;  // unit rotation axis
    double theta = std::atan2(s, c);
    double ct = std::cos(theta), st = std::sin(theta);
    // Rodrigues: R = I cosθ + sinθ [k]_x + (1-cosθ) k k^T
    Mat3 r;
    double K[3][3] = {{0, -k.z, k.y}, {k.z, 0, -k.x}, {-k.y, k.x, 0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r.m[i][j] = (i == j ? ct : 0.0) + st * K[i][j] + (1.0 - ct) * k[i] * k[j];
    return r;
}

}  // namespace msn
