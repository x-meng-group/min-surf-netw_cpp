#include "network.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <stdexcept>

#include "adgeom.hpp"
#include "autodiff.hpp"
#include "lbfgs.hpp"

namespace msn {

static const double PI = 3.14159265358979323846;
static int mround(double x) { return static_cast<int>(std::nearbyint(x)); }

double Network::regularPolygonPerimeter(int n) const {
    return n * 2.0 * std::sin(PI / n);
}

Vec3 Network::meanVertexPolyCoords(int v1) const {
    const auto& c = P.vertexPolyCoordinates[v1 - 1];
    Vec3 s{0, 0, 0};
    for (const auto& p : c) s += p;
    return s / static_cast<double>(c.size());
}

// face corner coordinate indices (1-based) for vertex v0, face faceN0 (0-based).
std::vector<int> Network::faceCoordIndices(int v0, int faceN0) const {
    std::vector<int> out;
    for (int s : P.vertexPolyFacesByLines[v0][faceN0]) {
        int m = std::abs(s);
        const auto& line = P.vertexPolyLines[v0][m - 1];  // {coordA, coordB} 1-based
        out.push_back(s > 0 ? line[0] : line[1]);         // Sign>0 -> first component, else last
    }
    return out;
}

std::vector<int> Network::computeFaceSideLengths(int v0, int faceN0) const {
    int edgeNSigned = P.vertexPolyConnectionsByEdgeN[v0][faceN0];
    const auto& face = P.vertexPolyFacesByLines[v0][faceN0];
    std::vector<int> lens;
    for (int s : face) {
        int m = std::abs(s);
        lens.push_back(mround(P.vertexPolyLineLengths[v0][m - 1] / dW_));
    }
    int W = width_[std::abs(edgeNSigned) - 1];
    int k = static_cast<int>(lens.size()) - 1;
    while (k >= 0 && lens[k] == 0) --k;
    int total = 0;
    for (int v : lens) total += v;
    if (k >= 0) lens[k] += W - total;
    return lens;
}

std::pair<int, int> Network::firstPositionConnByEdge(int signedEdge) const {
    for (int v = 0; v < P.nVertices; ++v)
        for (int f = 0; f < static_cast<int>(P.vertexPolyConnectionsByEdgeN[v].size()); ++f)
            if (P.vertexPolyConnectionsByEdgeN[v][f] == signedEdge) return {v, f};
    throw std::runtime_error("connection not found for signed edge");
}

// ---------------------------------------------------------------------------
// setupFEM
// ---------------------------------------------------------------------------
void Network::setupFEM(double globalScale) {
    const int E = P.numEdges();
    dW_ = 0.12566370614359174 * globalScale;
    dL_ = dW_;
    width_.assign(E, 0);
    length_.assign(E, 0);
    radius_.assign(E, 0.0);
    globalOffset_.assign(E, 0);
    quadOffset_.assign(E, 0);
    quadsLocal_.assign(E, {});
    boundaryMeshLabelsS_.assign(E, {});
    boundaryMeshLabelsT_.assign(E, {});
    numNodes_ = 0;
    numQuads_ = 0;

    for (int e = 0; e < E; ++e) {
        int w = mround(P.edgeWidthList[e] / dW_);
        if (w < 3) w = 3;
        int len = mround(P.edgeLengthList[e] / dL_);
        if (len < 2) len = 2;
        width_[e] = w;
        length_[e] = len;
        radius_[e] = P.edgeWidthList[e] / regularPolygonPerimeter(w);
        globalOffset_[e] = numNodes_;
        quadOffset_[e] = numQuads_;
        numNodes_ += w * len;
        numQuads_ += (len - 1) * w;

        auto nl = [&](int l, int ww) { return (ww - 1) * len + (l - 1); };  // local 0-based
        for (int ww = 1; ww <= w; ++ww) {
            int w2 = (ww % w) + 1;
            for (int l = 1; l <= len - 1; ++l)
                quadsLocal_[e].push_back({nl(l, ww), nl(l + 1, ww), nl(l + 1, w2), nl(l, w2)});
        }
        boundaryMeshLabelsS_[e].resize(w);
        boundaryMeshLabelsT_[e].resize(w);
        for (int ww = 1; ww <= w; ++ww) {
            boundaryMeshLabelsS_[e][ww - 1] = globalOffset_[e] + (ww - 1) * len + 0;
            boundaryMeshLabelsT_[e][ww - 1] = globalOffset_[e] + (ww - 1) * len + (len - 1);
        }
    }

    label_.assign(numNodes_, {});
    for (int e = 0; e < E; ++e) {
        int w = width_[e], len = length_[e];
        for (int ww = 1; ww <= w; ++ww) {
            double ang = 2.0 * PI * (ww - 1) / w;
            double rx = radius_[e] * std::cos(ang);
            double ry = radius_[e] * std::sin(ang);
            for (int l = 1; l <= len; ++l) {
                int g = globalOffset_[e] + (ww - 1) * len + (l - 1);
                label_[g] = {e, l, ww, rx, ry};
            }
        }
    }
    pbcOffset_.assign(numNodes_, Vec3{0, 0, 0});
}

// ---------------------------------------------------------------------------
// doAttribute
// ---------------------------------------------------------------------------
std::vector<int> Network::findBoundaryMeshLabelsPerFaceSide(int v0, int faceN0,
                                                            int faceSideNSigned) const {
    int edgeNSigned = P.vertexPolyConnectionsByEdgeN[v0][faceN0];
    if (edgeNSigned == 0) return {};
    std::vector<int> bml;
    int twist = 0;
    if (edgeNSigned > 0) {
        bml = boundaryMeshLabelsS_[edgeNSigned - 1];
    } else {
        int e = -edgeNSigned - 1;
        bml = boundaryMeshLabelsT_[e];
        std::rotate(bml.begin(), bml.begin() + 1, bml.end());  // RotateLeft 1
        std::reverse(bml.begin(), bml.end());                  // Reverse
        twist = mround(width_[e] * P.edgeTwistList[e] / (2 * PI));
    }
    const std::vector<int>& sideLens = vpFaceSideLengths_[v0][faceN0];

    // findLabelPartitions
    std::vector<int> acc(sideLens.size() + 1, 0);
    for (std::size_t i = 0; i < sideLens.size(); ++i) acc[i + 1] = acc[i] + sideLens[i];
    int total = acc.back();
    std::vector<int> range(3 * total);
    for (int i = 0; i < total; ++i) range[i] = range[i + total] = range[i + 2 * total] = i + 1;
    // RotateLeft by twist (handle negative)
    int t = ((twist % (3 * total)) + (3 * total)) % (3 * total);
    std::rotate(range.begin(), range.begin() + t, range.end());

    auto partitionR = [&](int r) {  // 0-based r
        std::vector<int> part;
        for (int i = total + acc[r]; i <= total + acc[r + 1]; ++i) part.push_back(range[i]);
        return part;
    };

    std::vector<int> sel;
    if (faceSideNSigned > 0) {
        sel = partitionR(faceSideNSigned - 1);
    } else {
        sel = partitionR(-faceSideNSigned - 1);
        std::reverse(sel.begin(), sel.end());
    }
    std::vector<int> out;
    for (int idx1 : sel) out.push_back(bml[idx1 - 1]);  // idx1 is 1-based into bml
    return out;
}

std::vector<Vec3> Network::generateBoundaryCoordinates(int signedEdge,
                                                       std::vector<int>& sideLensOut) const {
    auto pos = firstPositionConnByEdge(signedEdge);
    int v0 = pos.first, faceN0 = pos.second;
    std::vector<int> coordIdx = vpFacesByCoords_[v0][faceN0];
    std::vector<Vec3> polyCoordinates;
    for (int ci : coordIdx) polyCoordinates.push_back(P.vertexPolyCoordinates[v0][ci - 1]);
    std::vector<int> faceSideLengths = vpFaceSideLengths_[v0][faceN0];
    const std::vector<int> originalFace = P.vertexPolyFacesByLines[v0][faceN0];

    // minFaceSideLengths: min side length among all faces sharing each line.
    std::vector<int> minSide(originalFace.size());
    for (std::size_t e = 0; e < originalFace.size(); ++e) {
        int m = std::abs(originalFace[e]);
        int mn = faceSideLengths[e];
        for (std::size_t f2 = 0; f2 < P.vertexPolyFacesByLines[v0].size(); ++f2)
            for (std::size_t s2 = 0; s2 < P.vertexPolyFacesByLines[v0][f2].size(); ++s2)
                if (std::abs(P.vertexPolyFacesByLines[v0][f2][s2]) == m)
                    mn = std::min(mn, vpFaceSideLengths_[v0][f2][s2]);
        minSide[e] = mn;
    }

    if (signedEdge < 0) {
        std::rotate(polyCoordinates.begin(), polyCoordinates.begin() + 1, polyCoordinates.end());
        std::reverse(polyCoordinates.begin(), polyCoordinates.end());
        std::reverse(faceSideLengths.begin(), faceSideLengths.end());
        std::reverse(minSide.begin(), minSide.end());
    }
    std::vector<Vec3> polyEdgeVectors(polyCoordinates.size());
    for (std::size_t i = 0; i < polyCoordinates.size(); ++i)
        polyEdgeVectors[i] = polyCoordinates[(i + 1) % polyCoordinates.size()] - polyCoordinates[i];

    std::vector<Vec3> ring;
    sideLensOut.clear();
    for (std::size_t e = 0; e < polyEdgeVectors.size(); ++e) {
        std::vector<Vec3> base;
        if (minSide[e] == 0)
            base = {polyCoordinates[e]};
        else
            base = resampleSegment(polyCoordinates[e], polyCoordinates[e] + polyEdgeVectors[e],
                                   minSide[e] + 1);
        int pad = faceSideLengths[e] - minSide[e];
        std::vector<Vec3> padded;
        if (originalFace[e] > 0) {
            padded = base;
            for (int k = 0; k < pad; ++k) padded.push_back(base.back());
        } else {
            for (int k = 0; k < pad; ++k) padded.push_back(base.front());
            padded.insert(padded.end(), base.begin(), base.end());
        }
        padded.pop_back();  // Most
        for (const auto& p : padded) ring.push_back(p);
        sideLensOut.push_back(static_cast<int>(padded.size()));
    }
    return ring;
}

void Network::doAttribute() {
    const int E = P.numEdges();
    vpNumFaces_.assign(P.nVertices, 0);
    for (int v = 0; v < P.nVertices; ++v)
        vpNumFaces_[v] = static_cast<int>(P.vertexPolyConnectionsByEdgeN[v].size());

    edgeDisplacement_.assign(E, {});
    for (int e = 0; e < E; ++e)
        edgeDisplacement_[e] = P.edgePathList[e].back() - P.edgePathList[e].front();

    vpFacesByCoords_.assign(P.nVertices, {});
    vpFaceSideLengths_.assign(P.nVertices, {});
    for (int v = 0; v < P.nVertices; ++v) {
        int nf = static_cast<int>(P.vertexPolyFacesByLines[v].size());
        vpFacesByCoords_[v].resize(nf);
        vpFaceSideLengths_[v].resize(nf);
        for (int f = 0; f < nf; ++f) {
            vpFacesByCoords_[v][f] = faceCoordIndices(v, f);
            vpFaceSideLengths_[v][f] = computeFaceSideLengths(v, f);
        }
    }

    // path interpolation with leftPad/rightPad
    pathInterp_.assign(E, {});
    for (int e = 0; e < E; ++e) {
        const auto& path = P.edgePathList[e];
        std::vector<double> intrinsic(path.size(), 0.0);
        for (std::size_t i = 1; i < path.size(); ++i)
            intrinsic[i] = intrinsic[i - 1] + norm(path[i] - path[i - 1]);
        double last = intrinsic.back();

        auto padFor = [&](int signedEdge, int vertex1Based) {
            auto pos = firstPositionConnByEdge(signedEdge);
            std::vector<int> ci = vpFacesByCoords_[pos.first][pos.second];
            Vec3 sum{0, 0, 0};
            for (int c : ci) sum += P.vertexPolyCoordinates[vertex1Based - 1][c - 1];
            Vec3 meanFace = sum / static_cast<double>(ci.size());
            Vec3 meanAll = meanVertexPolyCoords(vertex1Based);
            double r = 0.2 * P.edgeWidthList[e] + norm(meanFace - meanAll);
            return (length_[e] - 1) / last * static_cast<double>(mround(r));
        };
        double leftPad = padFor(e + 1, P.edgeList[e][0]);
        double rightPad = padFor(-(e + 1), P.edgeList[e][1]);

        std::vector<double> knots(path.size());
        for (std::size_t i = 0; i < path.size(); ++i)
            knots[i] = intrinsic[i] / last * (length_[e] - 1 + leftPad + rightPad) + 1 - leftPad;
        pathInterp_[e] = LinearInterp(knots, path);
    }

    // parallel-transport rotation frames
    pathRot_.assign(E, {});
    for (int e = 0; e < E; ++e) {
        const LinearInterp& ip = pathInterp_[e];
        int len = length_[e];
        std::vector<Mat3> frames;
        frames.push_back(Mat3::identity());
        Vec3 normalPre = normalize(ip(1.0) - ip(0.0));
        Mat3 rotate = Mat3::identity();
        int steps = mround((len - 2) / 0.01) + 1;
        for (int k = 0; k < steps; ++k) {
            double tt = 2.0 + k * 0.01;
            Vec3 normalPost = normalize(ip(tt) - ip(tt - 1.0));
            if (vectorAngle(normalPre, normalPost) > 1e-10)
                rotate = rotationMatrixBetween(normalPre, normalPost) * rotate;
            normalPre = normalPost;
            if (k % 100 == 0) frames.push_back(rotate);
        }
        while (static_cast<int>(frames.size()) < len) frames.push_back(rotate);
        pathRot_[e] = frames;
    }

    // boundary rings + normals + fatter
    edgeBndS_.assign(E, {});
    edgeBndT_.assign(E, {});
    edgeBndFatS_.assign(E, {});
    edgeBndFatT_.assign(E, {});
    edgeBndNormalS_.assign(E, {});
    edgeBndNormalT_.assign(E, {});
    edgeBndSideLenS_.assign(E, {});
    edgeBndSideLenT_.assign(E, {});
    for (int e = 0; e < E; ++e) {
        edgeBndS_[e] = generateBoundaryCoordinates(e + 1, edgeBndSideLenS_[e]);
        edgeBndT_[e] = generateBoundaryCoordinates(-(e + 1), edgeBndSideLenT_[e]);

        auto computeNormal = [&](int vertex1Based, const std::vector<Vec3>& ring) {
            if (vpNumFaces_[vertex1Based - 1] <= 3) return normalize(edgeDisplacement_[e]);
            Vec3 m{0, 0, 0};
            for (const auto& p : ring) m += p;
            m = m / static_cast<double>(ring.size());
            return normalize(cross(ring[0] - m, ring[1] - m));
        };
        edgeBndNormalS_[e] = computeNormal(P.edgeList[e][0], edgeBndS_[e]);
        edgeBndNormalT_[e] = computeNormal(P.edgeList[e][1], edgeBndT_[e]);

        auto fatten = [&](const std::vector<Vec3>& ring, const Vec3& normal) {
            std::vector<Vec3> out(ring.size());
            int n = static_cast<int>(ring.size());
            for (int i = 0; i < n; ++i) {
                Vec3 d = ring[i] - ring[(i - 1 + n) % n];  // p - RotateRight(p,1)
                out[i] = ring[i] + 0.1 * P.edgeWidthList[e] * normalize(cross(d, normal));
            }
            return out;
        };
        edgeBndFatS_[e] = fatten(edgeBndS_[e], edgeBndNormalS_[e]);
        edgeBndFatT_[e] = fatten(edgeBndT_[e], edgeBndNormalT_[e]);
    }

    // PBC offsets (zero when Mean(T poly coords) == last path point)
    for (int e = 0; e < E; ++e) {
        Vec3 off = meanVertexPolyCoords(P.edgeList[e][1]) - P.edgePathList[e].back();
        for (int lab : boundaryMeshLabelsT_[e]) {
            pbcOffset_[lab] = off;
            pbcOffset_[lab - 1] = off;
        }
    }
}

// ---------------------------------------------------------------------------
// doConstraint  (assemble structural data for the energy terms)
// ---------------------------------------------------------------------------
void Network::doConstraint() {
    const int E = P.numEdges();

    // quads (global labels)
    quads_.clear();
    for (int e = 0; e < E; ++e) {
        int off = globalOffset_[e];
        int qg = quadOffset_[e];
        for (const auto& q : quadsLocal_[e]) {
            Quad Q;
            for (int i = 0; i < 4; ++i) Q.f[i] = off + q[i];
            Q.edge = e;
            Q.globalQuad = qg++;
            quads_.push_back(Q);
        }
    }

    // interior fairness triples
    fairTriples_.clear();
    for (int e = 0; e < E; ++e) {
        int len = length_[e], w = width_[e], off = globalOffset_[e];
        auto loc = [&](int l, int ww) { return off + (ww - 1) * len + (l - 1); };
        // along tube axis (each w strip): consecutive triples
        for (int ww = 1; ww <= w; ++ww)
            for (int l = 1; l <= len - 2; ++l)
                fairTriples_.push_back({loc(l, ww), loc(l + 1, ww), loc(l + 2, ww)});
        // around circumference (each interior l): cyclic triples, drop first/last l
        for (int l = 2; l <= len - 1; ++l)
            for (int ww = 1; ww <= w; ++ww) {
                int w1 = ww, w2 = (ww % w) + 1, w3 = ((ww + 1) % w) + 1;
                fairTriples_.push_back({loc(l, w1), loc(l, w2), loc(l, w3)});
            }
    }

    // gluing pairs
    gluePairs_.clear();
    for (int v = 0; v < P.nVertices; ++v) {
        int nLines = static_cast<int>(P.vertexPolyLines[v].size());
        for (int polyLineN = 1; polyLineN <= nLines; ++polyLineN) {
            std::vector<int> faceNList, faceSideNSignedList;
            for (int f = 0; f < static_cast<int>(P.vertexPolyFacesByLines[v].size()); ++f)
                for (int s = 0; s < static_cast<int>(P.vertexPolyFacesByLines[v][f].size()); ++s) {
                    if (P.vertexPolyFacesByLines[v][f][s] == polyLineN) {
                        faceNList.push_back(f);
                        faceSideNSignedList.push_back(s + 1);
                    }
                }
            for (int f = 0; f < static_cast<int>(P.vertexPolyFacesByLines[v].size()); ++f)
                for (int s = 0; s < static_cast<int>(P.vertexPolyFacesByLines[v][f].size()); ++s) {
                    if (P.vertexPolyFacesByLines[v][f][s] == -polyLineN) {
                        faceNList.push_back(f);
                        faceSideNSignedList.push_back(-(s + 1));
                    }
                }
            if (faceNList.size() < 2) continue;
            std::vector<int> l1 = findBoundaryMeshLabelsPerFaceSide(v, faceNList[0], faceSideNSignedList[0]);
            std::vector<int> l2 = findBoundaryMeshLabelsPerFaceSide(v, faceNList[1], faceSideNSignedList[1]);
            int m = std::min(l1.size(), l2.size());
            for (int i = 0; i < m; ++i) gluePairs_.push_back({l1[i], l2[i]});
        }
    }

    // boundary fairness triples (across seams)
    std::set<int> anyBoundaryS;
    for (int e = 0; e < E; ++e)
        for (int lab : boundaryMeshLabelsS_[e]) anyBoundaryS.insert(lab);
    fairbTriples_.clear();
    for (const auto& g : gluePairs_) {
        int a = g[0], b = g[1];
        int a0 = anyBoundaryS.count(a) ? a + 1 : a - 1;
        int b0 = anyBoundaryS.count(b) ? b + 1 : b - 1;
        fairbTriples_.push_back({a0, a, b0});
    }

    // terminals: for each terminal vertex, each edge-end touching it
    terminals_.clear();
    for (int v = 0; v < P.nVertices; ++v) {
        if (!P.vertexPolyIsTerminals[v]) continue;
        int vertex1 = v + 1;
        for (int e = 0; e < E; ++e) {
            for (int comp = 0; comp < 2; ++comp) {
                if (P.edgeList[e][comp] == vertex1) {
                    TerminalConn tc;
                    tc.labels = (comp == 0) ? boundaryMeshLabelsS_[e] : boundaryMeshLabelsT_[e];
                    tc.center = meanVertexPolyCoords(vertex1);
                    tc.edge = e;
                    terminals_.push_back(tc);
                }
            }
        }
    }

    // path-mean groups
    pathMeans_.clear();
    for (int e = 0; e < E; ++e) {
        int len = length_[e], w = width_[e], off = globalOffset_[e];
        for (int t = 1; t <= len; ++t) {
            PathMeanGroup g;
            g.edge = e;
            for (int ww = 1; ww <= w; ++ww) g.labels.push_back(off + (ww - 1) * len + (t - 1));
            g.target = pathInterp_[e](static_cast<double>(t));
            pathMeans_.push_back(g);
        }
    }
}

// ---------------------------------------------------------------------------
// doInitialize
// ---------------------------------------------------------------------------
void Network::doInitialize(double dLambdaInit, double sqrtZScaleInit) {
    const int E = P.numEdges();
    dlBase_ = 3 * numNodes_;
    zsBase_ = dlBase_ + numQuads_;
    params_.assign(zsBase_ + E, 0.0);

    for (int q = 0; q < numQuads_; ++q) params_[dlBase_ + q] = dLambdaInit;
    for (int e = 0; e < E; ++e) params_[zsBase_ + e] = sqrtZScaleInit;

    for (int e = 0; e < E; ++e) {
        int len = length_[e], w = width_[e];
        const LinearInterp& ip = pathInterp_[e];

        std::vector<Vec3> coordsS = edgeBndFatS_[e];
        std::vector<Vec3> coordsT = edgeBndFatT_[e];
        int rot = mround(width_[e] * P.edgeTwistList[e] / (2 * PI));
        rot = ((rot % w) + w) % w;
        std::rotate(coordsT.begin(), coordsT.begin() + rot, coordsT.end());

        Vec3 normalPre = normalize(ip(1.0) - ip(0.0));
        Vec3 normalPost = normalize(ip(static_cast<double>(len)) - ip(static_cast<double>(len - 1)));

        std::vector<Vec3> accum(len + 1);  // 1-based
        accum[1] = ip(1.0);
        for (int t = 2; t <= len; ++t)
            accum[t] = accum[t - 1] + norm(ip(t) - ip(t - 1.0)) * normalPre;

        auto meanOf = [](const std::vector<Vec3>& v) {
            Vec3 s{0, 0, 0};
            for (const auto& p : v) s += p;
            return s / static_cast<double>(v.size());
        };

        Mat3 rotateS = (vectorAngle(edgeBndNormalS_[e], normalPre) <= 1e-10)
                           ? Mat3::identity()
                           : rotationMatrixBetween(edgeBndNormalS_[e], normalPre);
        Vec3 meanS = meanOf(coordsS);
        std::vector<Vec3> SFake(coordsS.size());
        for (std::size_t i = 0; i < coordsS.size(); ++i)
            SFake[i] = rotateS * (coordsS[i] - meanS) + accum[1];

        Mat3 rotateT = (vectorAngle(edgeBndNormalT_[e], normalPost) <= 1e-10)
                           ? Mat3::identity()
                           : rotationMatrixBetween(edgeBndNormalT_[e], normalPost);
        Mat3 rotateT2 = pathRot_[e][len - 1].transpose() * rotateT;
        Vec3 meanT = meanOf(coordsT);
        std::vector<Vec3> TFake(coordsT.size());
        for (std::size_t i = 0; i < coordsT.size(); ++i)
            TFake[i] = rotateT2 * (coordsT[i] - meanT) + accum[len];

        // per-row rings, then place along the path with the transport frames
        for (int t = 1; t <= len; ++t) {
            double a = (len > 1) ? static_cast<double>(t - 1) / (len - 1) : 0.0;
            std::vector<Vec3> ring(w);
            for (int i = 0; i < w; ++i) ring[i] = SFake[i] + (TFake[i] - SFake[i]) * a;
            Vec3 meanRing = meanOf(ring);
            Vec3 pathPt = ip(static_cast<double>(t));
            for (int ww = 1; ww <= w; ++ww) {
                Vec3 c = pathRot_[e][t - 1] * (ring[ww - 1] - meanRing) + pathPt;
                int g = globalOffset_[e] + (ww - 1) * len + (t - 1);
                params_[3 * g + 0] = c.x;
                params_[3 * g + 1] = c.y;
                params_[3 * g + 2] = c.z;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// doFindMinimum
// ---------------------------------------------------------------------------
EnergyBreakdown Network::doFindMinimum(const SolveWeights& W) {
    const int E = P.numEdges();
    double sumLW = 0.0;
    for (int e = 0; e < E; ++e) sumLW += P.edgeLengthList[e] * P.edgeWidthList[e];
    double gFactor = 1e3 / sumLW;
    double totalEdgeLength = 0.0;
    for (int e = 0; e < E; ++e) totalEdgeLength += P.edgeLengthList[e];

    Tape tape;
    // parameter leaves in fixed order matching params_
    std::vector<AD> adX(numNodes_), adY(numNodes_), adZ(numNodes_);
    for (int L = 0; L < numNodes_; ++L) {
        adX[L] = tape.variable(params_[3 * L + 0]);
        adY[L] = tape.variable(params_[3 * L + 1]);
        adZ[L] = tape.variable(params_[3 * L + 2]);
    }
    std::vector<AD> adDl(numQuads_);
    for (int q = 0; q < numQuads_; ++q) adDl[q] = tape.variable(params_[dlBase_ + q]);
    std::vector<AD> adZs(E);
    for (int e = 0; e < E; ++e) adZs[e] = tape.variable(params_[zsBase_ + e]);

    auto coord = [&](int L) { return Vec3AD{adX[L], adY[L], adZ[L]}; };
    auto coordPBC = [&](int L) {
        const Vec3& o = pbcOffset_[L];
        return Vec3AD{adX[L] + o.x, adY[L] + o.y, adZ[L] + o.z};
    };
    // rest (intrinsic) coordinate; z depends on sqrtZScale of the node's edge
    auto restZ = [&](int L) {
        const auto& li = label_[L];
        return adsquare(adZs[li.edge]) * (dL_ * li.l);
    };
    auto restXYsq = [&](int L1, int L2) {
        double dx = label_[L1].restX - label_[L2].restX;
        double dy = label_[L1].restY - label_[L2].restY;
        return dx * dx + dy * dy;
    };

    AD energy = tape.constant(0.0);
    EnergyBreakdown bd;

    // 1a. isometry + 1b. conformal factors (and surface uses mean dλ²)
    std::vector<AD> sumDlSq(E);
    std::vector<int> quadCountPerEdge(E, 0);
    for (int e = 0; e < E; ++e) sumDlSq[e] = tape.constant(0.0);
    AD isoEnergy = tape.constant(0.0);
    for (const auto& Q : quads_) {
        int f1 = Q.f[0], f2 = Q.f[1], f3 = Q.f[2], f4 = Q.f[3];
        AD dlSq = adsquare(adDl[Q.globalQuad]);
        sumDlSq[Q.edge] = sumDlSq[Q.edge] + dlSq;
        quadCountPerEdge[Q.edge]++;

        // rest diagonals (13) and (24)
        AD rz13 = restZ(f1) - restZ(f3);
        AD rz24 = restZ(f2) - restZ(f4);
        AD restLen13 = adsquare(rz13) + restXYsq(f1, f3);
        AD restLen24 = adsquare(rz24) + restXYsq(f2, f4);
        double rxy13x = label_[f1].restX - label_[f3].restX;
        double rxy13y = label_[f1].restY - label_[f3].restY;
        double rxy24x = label_[f2].restX - label_[f4].restX;
        double rxy24y = label_[f2].restY - label_[f4].restY;
        AD restDot = rz13 * rz24 + (rxy13x * rxy24x + rxy13y * rxy24y);

        Vec3AD d13 = coord(f1) - coord(f3);
        Vec3AD d24 = coord(f2) - coord(f4);
        AD A = (dlSq + 1.0) * restLen13 - normSq(d13);
        AD B = (dlSq + 1.0) * restLen24 - normSq(d24);
        AD C = (dlSq + 1.0) * restDot - dot(d13, d24);
        isoEnergy = isoEnergy + (adsquare(A) + adsquare(B) + adsquare(C));
    }
    energy = energy + isoEnergy * W.wIso;

    // 1c. gluing
    AD glueEnergy = tape.constant(0.0);
    for (const auto& g : gluePairs_) {
        Vec3AD d = coordPBC(g[0]) - coordPBC(g[1]);
        glueEnergy = glueEnergy + normSq(d);
    }
    energy = energy + glueEnergy * 1e3;

    // 1d. interior fairness (scale-free)
    AD fairEnergy = tape.constant(0.0);
    for (const auto& tr : fairTriples_) {
        Vec3AD p1 = coord(tr[0]), p2 = coord(tr[1]), p3 = coord(tr[2]);
        Vec3AD sd = p1 - p2 * 2.0 + p3;
        AD denom = normSq(p1 - p2) + normSq(p2 - p3);
        AD inv = adpow(denom, -0.5);  // 1/sqrt(denom)
        fairEnergy = fairEnergy + (adsquare(sd.x) + adsquare(sd.y) + adsquare(sd.z)) * 1.0 *
                                      (inv * inv);
    }
    energy = energy + fairEnergy * W.wfair;

    // boundary fairness (scale-free)
    AD fairbEnergy = tape.constant(0.0);
    for (const auto& tr : fairbTriples_) {
        Vec3AD p1 = coordPBC(tr[0]), p2 = coordPBC(tr[1]), p3 = coordPBC(tr[2]);
        Vec3AD sd = p1 - p2 * 2.0 + p3;
        AD denom = normSq(p1 - p2) + normSq(p2 - p3);
        AD inv = adpow(denom, -0.5);
        fairbEnergy = fairbEnergy + (adsquare(sd.x) + adsquare(sd.y) + adsquare(sd.z)) * (inv * inv);
    }
    energy = energy + fairbEnergy * (1e-2 * totalEdgeLength * W.wfair);

    // 2c. terminal mean
    AD termMeanEnergy = tape.constant(0.0);
    for (const auto& tc : terminals_) {
        double w = static_cast<double>(width_[tc.edge]);
        Vec3AD sum{tape.constant(0.0), tape.constant(0.0), tape.constant(0.0)};
        for (int lab : tc.labels) sum = sum + coord(lab);
        AD rx = sum.x * (1.0 / w) - tc.center.x;
        AD ry = sum.y * (1.0 / w) - tc.center.y;
        AD rz = sum.z * (1.0 / w) - tc.center.z;
        termMeanEnergy = termMeanEnergy + adsquare(rx) + adsquare(ry) + adsquare(rz);
    }
    energy = energy + termMeanEnergy * 1e3;

    // 2d. terminal circle (push boundary away from center)
    AD termCircleEnergy = tape.constant(0.0);
    for (const auto& tc : terminals_) {
        double r2 = radius_[tc.edge] * radius_[tc.edge];
        AD inner = tape.constant(0.0);
        for (int lab : tc.labels) {
            Vec3AD d{adX[lab] - tc.center.x, adY[lab] - tc.center.y, adZ[lab] - tc.center.z};
            AD term = normSq(d) - r2;
            inner = inner + adsquare(term);
        }
        termCircleEnergy = termCircleEnergy + adsquare(inner);
    }
    energy = energy + termCircleEnergy * W.wCircle;

    // surface area term
    AD surfaceEnergy = tape.constant(0.0);
    for (int e = 0; e < E; ++e) {
        double coef = W.wSurface * P.edgeLengthList[e] * P.edgeWidthList[e];
        AD meanConf = (quadCountPerEdge[e] > 0)
                          ? (sumDlSq[e] * (1.0 / quadCountPerEdge[e]) + 1.0)
                          : tape.constant(1.0);
        surfaceEnergy = surfaceEnergy + meanConf * adsquare(adZs[e]) * coef;
    }
    energy = energy + surfaceEnergy;

    // 2b. path-mean (only if weighted)
    AD pathEnergy = tape.constant(0.0);
    if (W.wPathPerFair != 0.0) {
        for (const auto& g : pathMeans_) {
            double w = static_cast<double>(width_[g.edge]);
            Vec3AD sum{tape.constant(0.0), tape.constant(0.0), tape.constant(0.0)};
            for (int lab : g.labels) sum = sum + coord(lab);
            AD rx = sum.x * (1.0 / w) - g.target.x;
            AD ry = sum.y * (1.0 / w) - g.target.y;
            AD rz = sum.z * (1.0 / w) - g.target.z;
            pathEnergy = pathEnergy + adsquare(rx) + adsquare(ry) + adsquare(rz);
        }
        energy = energy + pathEnergy * (W.wPathPerFair * W.wfair);
    }

    AD total = energy * gFactor;

    // optimize
    std::vector<double> x = params_;
    Objective obj = [&](const std::vector<double>& p, std::vector<double>& grad) {
        tape.setParamsAndForward(p);
        return tape.backward(total.idx, grad);
    };
    LBFGSOptions opt;
    opt.maxIterations = W.maxIterations;
    opt.verbose = W.verbose;
    LBFGSResult res = lbfgs(obj, x, opt);
    params_ = x;

    // report energy breakdown at the solution
    tape.setParamsAndForward(params_);
    tape.forward();
    auto val = [&](AD a) { return tape.val[a.idx]; };
    bd.iso = gFactor * W.wIso * val(isoEnergy);
    bd.glue = gFactor * 1e3 * val(glueEnergy);
    bd.fairness = gFactor * W.wfair * val(fairEnergy);
    bd.boundaryFairness = gFactor * (1e-2 * totalEdgeLength * W.wfair) * val(fairbEnergy);
    bd.terminalMean = gFactor * 1e3 * val(termMeanEnergy);
    bd.terminalCircle = gFactor * W.wCircle * val(termCircleEnergy);
    bd.surface = gFactor * val(surfaceEnergy);
    bd.path = gFactor * (W.wPathPerFair * W.wfair) * val(pathEnergy);
    bd.total = val(total);
    if (W.verbose)
        std::printf("  energy total=%.8g iso=%.4g glue=%.4g fair=%.4g bfair=%.4g termMean=%.4g termCircle=%.4g surf=%.4g (lbfgs iters=%d |g|=%.3g)\n",
                    bd.total, bd.iso, bd.glue, bd.fairness, bd.boundaryFairness, bd.terminalMean,
                    bd.terminalCircle, bd.surface, res.iterations, res.gradNorm);
    return bd;
}

std::vector<std::vector<Vec3>> Network::edgeCoordinates() const {
    std::vector<std::vector<Vec3>> out(P.numEdges());
    for (int e = 0; e < P.numEdges(); ++e) {
        int n = width_[e] * length_[e];
        out[e].resize(n);
        for (int i = 0; i < n; ++i) {
            int g = globalOffset_[e] + i;
            out[e][i] = {params_[3 * g + 0], params_[3 * g + 1], params_[3 * g + 2]};
        }
    }
    return out;
}

}  // namespace msn
