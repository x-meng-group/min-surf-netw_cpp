// Output helpers: write the optimized network as Wavefront .obj meshes (one
// quad-mesh tube per edge) plus a small text dump of the per-edge parameters.
#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "network.hpp"

namespace msn {

// Write one .obj file containing every edge's quad-mesh tube. Quads are emitted
// as faces with 1-based vertex indices accumulated across edges.
inline void writeObj(const Network& net, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output file: " + path);
    f << "# min-surf-netw-cpp output mesh\n";

    std::vector<std::vector<Vec3>> coords = net.edgeCoordinates();
    const auto& quadsLocal = net.edgeQuadsLocal();

    int base = 0;  // running 0-based vertex count
    for (int e = 0; e < net.numEdges(); ++e) {
        f << "o edge_" << (e + 1) << "\n";
        for (const auto& v : coords[e])
            f << "v " << v.x << " " << v.y << " " << v.z << "\n";
        for (const auto& q : quadsLocal[e])
            f << "f " << (base + q[0] + 1) << " " << (base + q[1] + 1) << " "
              << (base + q[2] + 1) << " " << (base + q[3] + 1) << "\n";
        base += static_cast<int>(coords[e].size());
    }
}

// Dump per-edge solver parameters (sqrtZScale, mesh dimensions) to a text file.
inline void writeParamDump(const Network& net, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output file: " + path);
    f << "edge\twidth\tlength\tsqrtZScale\n";
    for (int e = 0; e < net.numEdges(); ++e)
        f << (e + 1) << "\t" << net.width(e) << "\t" << net.length(e) << "\t"
          << net.sqrtZScale(e) << "\n";
}

}  // namespace msn
