// Network: faithful C++ port of QuadMeshNetworkScript.wl.
//
// Pipeline (mirrors the Mathematica function names):
//   setupFEM      -> build a structured quad cylinder mesh per edge
//   doAttribute   -> path interpolation, parallel-transport frames, boundary rings
//   doConstraint  -> assemble the structural data for every energy term
//   doInitialize  -> place an initial 3D embedding along the guide paths
//   doFindMinimum -> minimize the energy with L-BFGS (replaces FindMinimum)
//
// The energy itself is built as a reverse-mode autodiff expression each time
// doFindMinimum runs, so the optimizer gets exact gradients.
#pragma once

#include <string>
#include <vector>

#include "interpolation.hpp"
#include "linalg.hpp"
#include "preset.hpp"

namespace msn {

struct SolveWeights {
    double wfair = 100.0;
    int maxIterations = 5000;
    double wSurface = 1e-3;
    double wCircle = 1e3;
    double wPathPerFair = 0.0;
    double wIso = 1.0;
    bool verbose = true;
};

struct EnergyBreakdown {
    double total = 0.0;
    double iso = 0.0;
    double fairness = 0.0;
    double boundaryFairness = 0.0;
    double surface = 0.0;
    double path = 0.0;
    double glue = 0.0;
    double terminalMean = 0.0;
    double terminalCircle = 0.0;
};

class Network {
public:
    explicit Network(const Preset& preset) : P(preset) {}

    void setupFEM(double globalScale = 1.0);
    void doAttribute();
    void doConstraint();
    void doInitialize(double dLambdaInit = 1.0, double sqrtZScaleInit = 1.0);
    // Runs one minimization at the given weights, warm-started from current state.
    EnergyBreakdown doFindMinimum(const SolveWeights& w);

    // Optimized coordinates per edge, ordered node(l,w) = (w-1)*length + (l-1).
    std::vector<std::vector<Vec3>> edgeCoordinates() const;
    // Quad connectivity per edge in terms of local node indices (0-based).
    const std::vector<std::vector<std::array<int, 4>>>& edgeQuadsLocal() const { return quadsLocal_; }

    int numEdges() const { return P.numEdges(); }
    int width(int e) const { return width_[e]; }
    int length(int e) const { return length_[e]; }
    double sqrtZScale(int e) const { return params_[zsBase_ + e]; }
    double dLambda(int globalQuad) const { return params_[dlBase_ + globalQuad]; }

private:
    // ---- preset ----
    Preset P;

    // ---- setupFEM outputs ----
    double dW_ = 0.0, dL_ = 0.0;
    std::vector<int> width_, length_;
    std::vector<double> radius_;
    std::vector<int> globalOffset_;     // first global label (0-based) of each edge
    std::vector<int> quadOffset_;       // first global quad index of each edge
    int numNodes_ = 0, numQuads_ = 0;
    std::vector<std::vector<std::array<int, 4>>> quadsLocal_;   // per edge, local node indices
    std::vector<std::array<int, 2>> boundaryS_, boundaryT_;     // [edge] -> not used directly; see below
    std::vector<std::vector<int>> boundaryMeshLabelsS_, boundaryMeshLabelsT_;  // global labels
    // per global label: edge, l(1..len), w(1..wid) and rest x,y
    struct LabelInfo { int edge, l, w; double restX, restY; };
    std::vector<LabelInfo> label_;

    // ---- doAttribute outputs ----
    std::vector<std::vector<int>> vpConnLen_;  // (unused placeholder)
    std::vector<int> vpNumFaces_;              // per vertex number of faces
    std::vector<Vec3> edgeDisplacement_;
    std::vector<std::vector<std::vector<int>>> vpFacesByCoords_;  // per vertex/face -> coord indices (1-based)
    std::vector<std::vector<std::vector<int>>> vpFaceSideLengths_;  // per vertex/face -> per-line side length
    std::vector<LinearInterp> pathInterp_;
    std::vector<std::vector<Mat3>> pathRot_;  // per edge, frames for l=1..length (index 0..length-1)
    // boundary rings per edge: per face-side list of Vec3 (flattened S/T) and fatter versions
    std::vector<std::vector<Vec3>> edgeBndS_, edgeBndT_, edgeBndFatS_, edgeBndFatT_;
    std::vector<Vec3> edgeBndNormalS_, edgeBndNormalT_;
    std::vector<std::vector<int>> edgeBndSideLenS_, edgeBndSideLenT_;  // length of each face-side group
    std::vector<Vec3> pbcOffset_;  // per global label (0 default)

    // ---- doConstraint outputs ----
    struct Quad { int f[4]; int edge; int globalQuad; };
    std::vector<Quad> quads_;                       // all quads (global labels)
    std::vector<std::array<int, 3>> fairTriples_;   // per-edge interior fairness (global labels)
    std::vector<std::array<int, 3>> fairbTriples_;  // boundary fairness across seams
    std::vector<std::array<int, 2>> gluePairs_;     // (label1,label2)
    struct TerminalConn { std::vector<int> labels; Vec3 center; int edge; };
    std::vector<TerminalConn> terminals_;
    struct PathMeanGroup { int edge; std::vector<int> labels; Vec3 target; };
    std::vector<PathMeanGroup> pathMeans_;

    // ---- parameters (persistent across solves) ----
    std::vector<double> params_;
    int dlBase_ = 0, zsBase_ = 0;

    // helpers
    double regularPolygonPerimeter(int n) const;
    Vec3 meanVertexPolyCoords(int v1) const;  // v1 is 1-based
    std::vector<int> faceCoordIndices(int v0, int faceN0) const;  // 0-based -> 1-based coord idx list
    std::vector<int> computeFaceSideLengths(int v0, int faceN0) const;
    std::pair<int, int> firstPositionConnByEdge(int signedEdge) const;  // returns {v0, faceN0}
    std::vector<Vec3> generateBoundaryCoordinates(int signedEdge, std::vector<int>& sideLensOut) const;
    std::vector<int> findBoundaryMeshLabelsPerFaceSide(int v0, int faceN0, int faceSideNSigned) const;
};

}  // namespace msn
