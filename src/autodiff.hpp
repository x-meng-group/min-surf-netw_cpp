// Reverse-mode automatic differentiation over a static expression tape.
//
// The energy functional in the original Mathematica code is a single symbolic
// expression that FindMinimum differentiates analytically. Here we build the
// expression graph once as a tape of elementary operations, then evaluate value
// and full gradient by a forward + reverse sweep. Because the graph topology is
// fixed across solver iterations, we build it a single time and only update the
// leaf (parameter) values each evaluation.
#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

namespace msn {

enum class Op : unsigned char {
    Const,   // constant leaf (value fixed at creation)
    Var,     // parameter leaf (value updated each evaluation)
    Add,     // a + b
    Sub,     // a - b
    Mul,     // a * b
    Div,     // a / b
    Neg,     // -a
    AddC,    // a + c        (c in data)
    MulC,    // a * c        (c in data)
    Sqrt,    // sqrt(a)
    Sin,     // sin(a)
    Cos,     // cos(a)
    Square,  // a^2
    PowC,    // a^c          (c in data)
};

class Tape;

// Lightweight handle to a node in a Tape. Supports natural operator syntax for
// building expressions.
struct AD {
    Tape* tape = nullptr;
    int idx = -1;

    AD() = default;
    AD(Tape* t, int i) : tape(t), idx(i) {}

    AD operator+(const AD& o) const;
    AD operator-(const AD& o) const;
    AD operator*(const AD& o) const;
    AD operator/(const AD& o) const;
    AD operator-() const;
    AD operator+(double c) const;
    AD operator-(double c) const;
    AD operator*(double c) const;
};

class Tape {
public:
    std::vector<Op> op;
    std::vector<int> a;
    std::vector<int> b;
    std::vector<double> data;  // constant payload for Const/AddC/MulC/PowC
    std::vector<double> val;
    std::vector<double> adj;
    std::vector<int> varNodes;  // node indices of parameter leaves, in order

    int newNode(Op o, int aa, int bb, double c) {
        op.push_back(o);
        a.push_back(aa);
        b.push_back(bb);
        data.push_back(c);
        val.push_back(0.0);
        adj.push_back(0.0);
        return static_cast<int>(op.size()) - 1;
    }

    AD constant(double c) { return AD(this, newNode(Op::Const, -1, -1, c)); }

    // Creates a parameter leaf with an initial value.
    AD variable(double initial) {
        int i = newNode(Op::Var, -1, -1, initial);
        varNodes.push_back(i);
        return AD(this, i);
    }

    std::size_t numVars() const { return varNodes.size(); }
    std::size_t numNodes() const { return op.size(); }

    // Update parameter leaf values from a flat vector (same order as creation).
    // Values live in `data` for Var nodes; forward() copies them into `val`.
    void setParams(const std::vector<double>& p) {
        for (std::size_t i = 0; i < varNodes.size(); ++i) data[varNodes[i]] = p[i];
    }

    void forward() {
        const std::size_t n = op.size();
        for (std::size_t i = 0; i < n; ++i) {
            switch (op[i]) {
                case Op::Const:
                case Op::Var: val[i] = data[i]; break;  // Var: data holds current value
                case Op::Add: val[i] = val[a[i]] + val[b[i]]; break;
                case Op::Sub: val[i] = val[a[i]] - val[b[i]]; break;
                case Op::Mul: val[i] = val[a[i]] * val[b[i]]; break;
                case Op::Div: val[i] = val[a[i]] / val[b[i]]; break;
                case Op::Neg: val[i] = -val[a[i]]; break;
                case Op::AddC: val[i] = val[a[i]] + data[i]; break;
                case Op::MulC: val[i] = val[a[i]] * data[i]; break;
                case Op::Sqrt: val[i] = std::sqrt(val[a[i]]); break;
                case Op::Sin: val[i] = std::sin(val[a[i]]); break;
                case Op::Cos: val[i] = std::cos(val[a[i]]); break;
                case Op::Square: { double v = val[a[i]]; val[i] = v * v; } break;
                case Op::PowC: val[i] = std::pow(val[a[i]], data[i]); break;
            }
        }
    }

    void setParamsAndForward(const std::vector<double>& p) {
        setParams(p);
        forward();
    }

    // Reverse sweep from `root`; fills `grad` (size numVars()).
    double backward(int root, std::vector<double>& grad) {
        std::fill(adj.begin(), adj.end(), 0.0);
        adj[root] = 1.0;
        for (int i = static_cast<int>(op.size()) - 1; i >= 0; --i) {
            double g = adj[i];
            if (g == 0.0) continue;
            switch (op[i]) {
                case Op::Const:
                case Op::Var: break;
                case Op::Add: adj[a[i]] += g; adj[b[i]] += g; break;
                case Op::Sub: adj[a[i]] += g; adj[b[i]] -= g; break;
                case Op::Mul:
                    adj[a[i]] += g * val[b[i]];
                    adj[b[i]] += g * val[a[i]];
                    break;
                case Op::Div: {
                    double vb = val[b[i]];
                    adj[a[i]] += g / vb;
                    adj[b[i]] -= g * val[a[i]] / (vb * vb);
                } break;
                case Op::Neg: adj[a[i]] -= g; break;
                case Op::AddC: adj[a[i]] += g; break;
                case Op::MulC: adj[a[i]] += g * data[i]; break;
                case Op::Sqrt: {
                    double v = val[i];
                    if (v > 1e-300) adj[a[i]] += g * 0.5 / v;
                } break;
                case Op::Sin: adj[a[i]] += g * std::cos(val[a[i]]); break;
                case Op::Cos: adj[a[i]] -= g * std::sin(val[a[i]]); break;
                case Op::Square: adj[a[i]] += g * 2.0 * val[a[i]]; break;
                case Op::PowC: {
                    double c = data[i];
                    adj[a[i]] += g * c * std::pow(val[a[i]], c - 1.0);
                } break;
            }
        }
        grad.assign(numVars(), 0.0);
        for (std::size_t k = 0; k < varNodes.size(); ++k) grad[k] = adj[varNodes[k]];
        return val[root];
    }
};

inline AD AD::operator+(const AD& o) const { return AD(tape, tape->newNode(Op::Add, idx, o.idx, 0)); }
inline AD AD::operator-(const AD& o) const { return AD(tape, tape->newNode(Op::Sub, idx, o.idx, 0)); }
inline AD AD::operator*(const AD& o) const { return AD(tape, tape->newNode(Op::Mul, idx, o.idx, 0)); }
inline AD AD::operator/(const AD& o) const { return AD(tape, tape->newNode(Op::Div, idx, o.idx, 0)); }
inline AD AD::operator-() const { return AD(tape, tape->newNode(Op::Neg, idx, -1, 0)); }
inline AD AD::operator+(double c) const { return AD(tape, tape->newNode(Op::AddC, idx, -1, c)); }
inline AD AD::operator-(double c) const { return AD(tape, tape->newNode(Op::AddC, idx, -1, -c)); }
inline AD AD::operator*(double c) const { return AD(tape, tape->newNode(Op::MulC, idx, -1, c)); }

inline AD operator*(double c, const AD& v) { return v * c; }
inline AD operator+(double c, const AD& v) { return v + c; }
inline AD operator-(double c, const AD& v) { return (-v) + c; }

inline AD adsqrt(const AD& v) { return AD(v.tape, v.tape->newNode(Op::Sqrt, v.idx, -1, 0)); }
inline AD adsin(const AD& v) { return AD(v.tape, v.tape->newNode(Op::Sin, v.idx, -1, 0)); }
inline AD adcos(const AD& v) { return AD(v.tape, v.tape->newNode(Op::Cos, v.idx, -1, 0)); }
inline AD adsquare(const AD& v) { return AD(v.tape, v.tape->newNode(Op::Square, v.idx, -1, 0)); }
inline AD adpow(const AD& v, double c) { return AD(v.tape, v.tape->newNode(Op::PowC, v.idx, -1, c)); }

}  // namespace msn
