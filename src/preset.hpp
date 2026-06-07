// Preset = the topological + geometric description of a network manifold.
//
// This mirrors the variables that the original Mathematica `preset[]` builds and
// dumps to a `.mx` file (edgeList, vertexPolyLines, vertexPolyCoordinates, ...).
// To keep the JSON a faithful transcription of the `.mx` data, indices follow the
// original 1-based, signed convention:
//   * edge / vertex / line / coordinate indices are 1-based;
//   * signed edge ids in `vertexPolyConnectionsByEdgeN`: + = this vertex is the
//     S (source) end of the edge, - = the T (target) end, 0 = no end;
//   * signed line ids in `vertexPolyFacesByLines`: sign encodes orientation.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "json.hpp"
#include "linalg.hpp"

namespace msn {

struct Preset {
    std::string name;
    int nVertices = 0;

    std::vector<std::array<int, 2>> edgeList;  // 1-based {sourceVertex, targetVertex}
    std::vector<double> edgeLengthList;
    std::vector<double> edgeWidthList;
    std::vector<double> edgeTwistList;
    std::vector<std::vector<Vec3>> edgePathList;  // per edge: guide-path points

    std::vector<std::vector<std::array<int, 2>>> vertexPolyLines;  // per vertex/line: {coordA,coordB} 1-based
    std::vector<std::vector<double>> vertexPolyLineLengths;
    std::vector<std::vector<std::vector<int>>> vertexPolyFacesByLines;  // per vertex/face: signed line ids
    std::vector<std::vector<int>> vertexPolyConnectionsByEdgeN;         // per vertex: signed edge ids
    std::vector<std::vector<Vec3>> vertexPolyCoordinates;              // per vertex: polygon corners
    std::vector<bool> vertexPolyIsTerminals;

    int numEdges() const { return static_cast<int>(edgeList.size()); }

    void validate() const;
};

// Parse a Preset from JSON text.
Preset presetFromJson(const Json& j);
Preset loadPresetFile(const std::string& path);

// Serialize a Preset to JSON text (used to emit the built-in examples).
std::string presetToJson(const Preset& p);

// Built-in canonical presets (analogues of the Mathematica preset[] cases).
//   bifurcation  : 1 junction + 3 terminals (4 vertices, 3 edges)
//   trifurcation : 2 junctions + 4 terminals (6 vertices, 5 edges)
Preset makeBifurcation(double widthBimodalRatio = 1.0);
Preset makeTrifurcation(double lengthScale = 5.0);

}  // namespace msn
