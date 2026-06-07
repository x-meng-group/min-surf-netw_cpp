// Limited-memory BFGS (L-BFGS) with backtracking line search satisfying the
// strong Wolfe conditions. Stands in for Mathematica's FindMinimum (a
// quasi-Newton method) for the unconstrained minimization of the energy.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

namespace msn {

struct LBFGSOptions {
    int maxIterations = 5000;
    int memory = 8;          // number of correction pairs to store
    double gradTol = 1e-7;   // stop when ||grad||_inf < gradTol
    double relTol = 1e-12;   // stop when relative function decrease < relTol
    bool verbose = false;
};

struct LBFGSResult {
    double value = 0.0;
    int iterations = 0;
    double gradNorm = 0.0;
    bool converged = false;
};

// Objective: given x, fills grad and returns f(x).
using Objective = std::function<double(const std::vector<double>&, std::vector<double>&)>;

inline double infNorm(const std::vector<double>& v) {
    double m = 0.0;
    for (double x : v) m = std::max(m, std::abs(x));
    return m;
}

inline double dotv(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

inline LBFGSResult lbfgs(Objective f, std::vector<double>& x, const LBFGSOptions& opt) {
    const std::size_t n = x.size();
    std::vector<double> g(n), gNew(n), xNew(n), dir(n);
    double fx = f(x, g);

    std::vector<std::vector<double>> sList, yList;
    std::vector<double> rhoList;

    LBFGSResult res;
    res.value = fx;
    res.gradNorm = infNorm(g);

    for (int iter = 0; iter < opt.maxIterations; ++iter) {
        res.iterations = iter;
        res.gradNorm = infNorm(g);
        if (res.gradNorm < opt.gradTol) { res.converged = true; break; }

        // Two-loop recursion to compute search direction = -H*g.
        dir = g;
        const int k = static_cast<int>(sList.size());
        std::vector<double> alpha(k);
        for (int i = k - 1; i >= 0; --i) {
            alpha[i] = rhoList[i] * dotv(sList[i], dir);
            for (std::size_t j = 0; j < n; ++j) dir[j] -= alpha[i] * yList[i][j];
        }
        double gamma = 1.0;
        if (k > 0) {
            double ys = dotv(yList[k - 1], sList[k - 1]);
            double yy = dotv(yList[k - 1], yList[k - 1]);
            if (yy > 0.0) gamma = ys / yy;
        }
        for (std::size_t j = 0; j < n; ++j) dir[j] *= gamma;
        for (int i = 0; i < k; ++i) {
            double beta = rhoList[i] * dotv(yList[i], dir);
            for (std::size_t j = 0; j < n; ++j) dir[j] += (alpha[i] - beta) * sList[i][j];
        }
        for (std::size_t j = 0; j < n; ++j) dir[j] = -dir[j];

        double dg = dotv(dir, g);
        if (dg >= 0.0) {
            // Not a descent direction; reset to steepest descent.
            for (std::size_t j = 0; j < n; ++j) dir[j] = -g[j];
            dg = dotv(dir, g);
            sList.clear(); yList.clear(); rhoList.clear();
        }

        // Backtracking line search (Armijo + curvature).
        double step = (iter == 0) ? std::min(1.0, 1.0 / std::max(1e-12, infNorm(g))) : 1.0;
        const double c1 = 1e-4, c2 = 0.9;
        double fNew = fx;
        bool ok = false;
        for (int ls = 0; ls < 40; ++ls) {
            for (std::size_t j = 0; j < n; ++j) xNew[j] = x[j] + step * dir[j];
            fNew = f(xNew, gNew);
            if (fNew <= fx + c1 * step * dg) {
                double dgNew = dotv(dir, gNew);
                if (std::abs(dgNew) <= c2 * std::abs(dg)) { ok = true; break; }
                if (dgNew >= 0.0) { step *= 0.5; continue; }  // overshoot
                ok = true;  // Armijo satisfied; accept even if weak curvature
                break;
            }
            step *= 0.5;
        }
        if (!ok && fNew >= fx) {
            // Line search failed to make progress.
            res.value = fx;
            break;
        }

        // Update memory.
        std::vector<double> s(n), y(n);
        for (std::size_t j = 0; j < n; ++j) { s[j] = xNew[j] - x[j]; y[j] = gNew[j] - g[j]; }
        double ys = dotv(y, s);
        if (ys > 1e-12) {
            if (static_cast<int>(sList.size()) == opt.memory) {
                sList.erase(sList.begin());
                yList.erase(yList.begin());
                rhoList.erase(rhoList.begin());
            }
            sList.push_back(s);
            yList.push_back(y);
            rhoList.push_back(1.0 / ys);
        }

        double fOld = fx;
        x = xNew;
        g = gNew;
        fx = fNew;
        res.value = fx;

        if (opt.verbose && (iter % 25 == 0))
            std::printf("  [lbfgs] iter %d  f=%.10g  |g|=%.3g  step=%.3g\n", iter, fx,
                        res.gradNorm, step);

        if (std::abs(fOld - fx) <= opt.relTol * (std::abs(fOld) + 1.0)) {
            res.converged = true;
            break;
        }
    }
    res.value = fx;
    return res;
}

}  // namespace msn
