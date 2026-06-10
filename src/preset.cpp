#include "preset.hpp"

#include <cmath>
#include <fstream>
#include <sstream>

namespace msn {

static const double PI = 3.14159265358979323846;

void Preset::validate() const {
    auto fail = [](const std::string& m) { throw std::runtime_error("Preset invalid: " + m); };
    if (nVertices <= 0) fail("nVertices <= 0");
    const int E = numEdges();
    if (static_cast<int>(edgeLengthList.size()) != E) fail("edgeLengthList size mismatch");
    if (static_cast<int>(edgeWidthList.size()) != E) fail("edgeWidthList size mismatch");
    if (static_cast<int>(edgeTwistList.size()) != E) fail("edgeTwistList size mismatch");
    if (static_cast<int>(edgePathList.size()) != E) fail("edgePathList size mismatch");
    if (static_cast<int>(vertexPolyLines.size()) != nVertices) fail("vertexPolyLines size");
    if (static_cast<int>(vertexPolyLineLengths.size()) != nVertices) fail("vertexPolyLineLengths size");
    if (static_cast<int>(vertexPolyFacesByLines.size()) != nVertices) fail("vertexPolyFacesByLines size");
    if (static_cast<int>(vertexPolyConnectionsByEdgeN.size()) != nVertices) fail("vertexPolyConnectionsByEdgeN size");
    if (static_cast<int>(vertexPolyCoordinates.size()) != nVertices) fail("vertexPolyCoordinates size");
    if (static_cast<int>(vertexPolyIsTerminals.size()) != nVertices) fail("vertexPolyIsTerminals size");
    for (const auto& e : edgeList)
        if (e[0] < 1 || e[0] > nVertices || e[1] < 1 || e[1] > nVertices) fail("edge vertex out of range");
}

// ----------------------------------------------------------------------------
// JSON (de)serialization
// ----------------------------------------------------------------------------
static Vec3 vec3FromJson(const Json& j) {
    return {j[0].asNumber(), j[1].asNumber(), j[2].asNumber()};
}

Preset presetFromJson(const Json& j) {
    Preset p;
    if (j.has("name")) p.name = j.at("name").asString();
    p.nVertices = j.at("nVertices").asInt();

    for (const auto& e : j.at("edgeList").arrayValue)
        p.edgeList.push_back({e[0].asInt(), e[1].asInt()});
    for (const auto& x : j.at("edgeLengthList").arrayValue) p.edgeLengthList.push_back(x.asNumber());
    for (const auto& x : j.at("edgeWidthList").arrayValue) p.edgeWidthList.push_back(x.asNumber());
    for (const auto& x : j.at("edgeTwistList").arrayValue) p.edgeTwistList.push_back(x.asNumber());

    for (const auto& path : j.at("edgePathList").arrayValue) {
        std::vector<Vec3> pts;
        for (const auto& pt : path.arrayValue) pts.push_back(vec3FromJson(pt));
        p.edgePathList.push_back(std::move(pts));
    }
    for (const auto& v : j.at("vertexPolyLines").arrayValue) {
        std::vector<std::array<int, 2>> lines;
        for (const auto& ln : v.arrayValue) lines.push_back({ln[0].asInt(), ln[1].asInt()});
        p.vertexPolyLines.push_back(std::move(lines));
    }
    for (const auto& v : j.at("vertexPolyLineLengths").arrayValue) {
        std::vector<double> lens;
        for (const auto& x : v.arrayValue) lens.push_back(x.asNumber());
        p.vertexPolyLineLengths.push_back(std::move(lens));
    }
    for (const auto& v : j.at("vertexPolyFacesByLines").arrayValue) {
        std::vector<std::vector<int>> faces;
        for (const auto& f : v.arrayValue) {
            std::vector<int> face;
            for (const auto& s : f.arrayValue) face.push_back(s.asInt());
            faces.push_back(std::move(face));
        }
        p.vertexPolyFacesByLines.push_back(std::move(faces));
    }
    for (const auto& v : j.at("vertexPolyConnectionsByEdgeN").arrayValue) {
        std::vector<int> conn;
        for (const auto& s : v.arrayValue) conn.push_back(s.asInt());
        p.vertexPolyConnectionsByEdgeN.push_back(std::move(conn));
    }
    for (const auto& v : j.at("vertexPolyCoordinates").arrayValue) {
        std::vector<Vec3> coords;
        for (const auto& c : v.arrayValue) coords.push_back(vec3FromJson(c));
        p.vertexPolyCoordinates.push_back(std::move(coords));
    }
    for (const auto& v : j.at("vertexPolyIsTerminals").arrayValue)
        p.vertexPolyIsTerminals.push_back(v.asBool());

    p.validate();
    return p;
}

Preset loadPresetFile(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open preset file: " + path);
    std::stringstream ss;
    ss << f.rdbuf();
    return presetFromJson(Json::parse(ss.str()));
}

static std::string num(double v) {
    std::ostringstream os;
    os.precision(17);
    os << v;
    return os.str();
}

std::string presetToJson(const Preset& p) {
    std::ostringstream o;
    o << "{\n";
    o << "  \"name\": \"" << p.name << "\",\n";
    o << "  \"nVertices\": " << p.nVertices << ",\n";

    o << "  \"edgeList\": [";
    for (std::size_t i = 0; i < p.edgeList.size(); ++i)
        o << (i ? ", " : "") << "[" << p.edgeList[i][0] << ", " << p.edgeList[i][1] << "]";
    o << "],\n";

    auto dblArr = [&](const char* key, const std::vector<double>& v, const char* end) {
        o << "  \"" << key << "\": [";
        for (std::size_t i = 0; i < v.size(); ++i) o << (i ? ", " : "") << num(v[i]);
        o << "]" << end;
    };
    dblArr("edgeLengthList", p.edgeLengthList, ",\n");
    dblArr("edgeWidthList", p.edgeWidthList, ",\n");
    dblArr("edgeTwistList", p.edgeTwistList, ",\n");

    o << "  \"edgePathList\": [\n";
    for (std::size_t e = 0; e < p.edgePathList.size(); ++e) {
        o << "    [";
        for (std::size_t i = 0; i < p.edgePathList[e].size(); ++i) {
            const Vec3& v = p.edgePathList[e][i];
            o << (i ? ", " : "") << "[" << num(v.x) << ", " << num(v.y) << ", " << num(v.z) << "]";
        }
        o << "]" << (e + 1 < p.edgePathList.size() ? ",\n" : "\n");
    }
    o << "  ],\n";

    o << "  \"vertexPolyLines\": [\n";
    for (std::size_t v = 0; v < p.vertexPolyLines.size(); ++v) {
        o << "    [";
        for (std::size_t i = 0; i < p.vertexPolyLines[v].size(); ++i)
            o << (i ? ", " : "") << "[" << p.vertexPolyLines[v][i][0] << ", "
              << p.vertexPolyLines[v][i][1] << "]";
        o << "]" << (v + 1 < p.vertexPolyLines.size() ? ",\n" : "\n");
    }
    o << "  ],\n";

    o << "  \"vertexPolyLineLengths\": [\n";
    for (std::size_t v = 0; v < p.vertexPolyLineLengths.size(); ++v) {
        o << "    [";
        for (std::size_t i = 0; i < p.vertexPolyLineLengths[v].size(); ++i)
            o << (i ? ", " : "") << num(p.vertexPolyLineLengths[v][i]);
        o << "]" << (v + 1 < p.vertexPolyLineLengths.size() ? ",\n" : "\n");
    }
    o << "  ],\n";

    o << "  \"vertexPolyFacesByLines\": [\n";
    for (std::size_t v = 0; v < p.vertexPolyFacesByLines.size(); ++v) {
        o << "    [";
        for (std::size_t f = 0; f < p.vertexPolyFacesByLines[v].size(); ++f) {
            o << (f ? ", " : "") << "[";
            for (std::size_t i = 0; i < p.vertexPolyFacesByLines[v][f].size(); ++i)
                o << (i ? ", " : "") << p.vertexPolyFacesByLines[v][f][i];
            o << "]";
        }
        o << "]" << (v + 1 < p.vertexPolyFacesByLines.size() ? ",\n" : "\n");
    }
    o << "  ],\n";

    o << "  \"vertexPolyConnectionsByEdgeN\": [";
    for (std::size_t v = 0; v < p.vertexPolyConnectionsByEdgeN.size(); ++v) {
        o << (v ? ", " : "") << "[";
        for (std::size_t i = 0; i < p.vertexPolyConnectionsByEdgeN[v].size(); ++i)
            o << (i ? ", " : "") << p.vertexPolyConnectionsByEdgeN[v][i];
        o << "]";
    }
    o << "],\n";

    o << "  \"vertexPolyCoordinates\": [\n";
    for (std::size_t v = 0; v < p.vertexPolyCoordinates.size(); ++v) {
        o << "    [";
        for (std::size_t i = 0; i < p.vertexPolyCoordinates[v].size(); ++i) {
            const Vec3& c = p.vertexPolyCoordinates[v][i];
            o << (i ? ", " : "") << "[" << num(c.x) << ", " << num(c.y) << ", " << num(c.z) << "]";
        }
        o << "]" << (v + 1 < p.vertexPolyCoordinates.size() ? ",\n" : "\n");
    }
    o << "  ],\n";

    o << "  \"vertexPolyIsTerminals\": [";
    for (std::size_t v = 0; v < p.vertexPolyIsTerminals.size(); ++v)
        o << (v ? ", " : "") << (p.vertexPolyIsTerminals[v] ? "true" : "false");
    o << "]\n";

    o << "}\n";
    return o.str();
}

// ----------------------------------------------------------------------------
// Built-in canonical presets
// ----------------------------------------------------------------------------
static Vec3 meanCoords(const std::vector<Vec3>& v) {
    Vec3 s{0, 0, 0};
    for (const auto& p : v) s += p;
    return s / static_cast<double>(v.size());
}

static void fillStraightPaths(Preset& p) {
    p.edgePathList.resize(p.numEdges());
    for (int e = 0; e < p.numEdges(); ++e) {
        const Vec3 a = meanCoords(p.vertexPolyCoordinates[p.edgeList[e][0] - 1]);
        const Vec3 b = meanCoords(p.vertexPolyCoordinates[p.edgeList[e][1] - 1]);
        p.edgePathList[e] = {a, b};
    }
}

// Bifurcation: 1 junction (vertex 1) + 3 terminals (vertices 2,3,4).
// Transcription of the Mathematica "3C03_bimodal" preset.
Preset makeBifurcation(double r) {
    Preset p;
    p.name = "bifurcation";
    p.nVertices = 4;
    p.edgeList = {{1, 2}, {1, 3}, {1, 4}};
    p.edgeLengthList = {3.0, 3.0, 3.0};  // 3/5 * {5,5,5}
    p.edgeWidthList = {2 * PI, 2 * PI, 2 * PI * r};
    p.edgeTwistList = {0.0, 0.0, 0.0};

    p.vertexPolyLines = {{{1, 2}, {1, 2}, {1, 2}}, {{1, 2}, {1, 2}}, {{1, 2}, {1, 2}}, {{1, 2}, {1, 2}}};
    p.vertexPolyLineLengths = {{2 * PI * r / 2, 2 * PI * (1 - r / 2), 2 * PI * r / 2},
                               {2 * PI * r / 2, 2 * PI * (1 - r / 2)},
                               {2 * PI * (1 - r / 2), 2 * PI * r / 2},
                               {2 * PI * r / 2, 2 * PI * r / 2}};
    p.vertexPolyFacesByLines = {{{1, -2}, {2, -3}, {3, -1}}, {{2, -1}}, {{2, -1}}, {{2, -1}}};
    p.vertexPolyConnectionsByEdgeN = {{1, 2, 3}, {-1}, {-2}, {-3}};

    const double L = p.edgeLengthList[0];
    const Vec3 zc{0, 0, PI / 2};
    auto pair = [](const Vec3& base, const Vec3& off) {
        return std::vector<Vec3>{base + off, base - off};
    };
    p.vertexPolyCoordinates = {
        pair(L * Vec3{1.0 / 4, 3.0 / 4, 0}, zc),
        pair(L * Vec3{1, 0, 0}, zc),
        pair(L * Vec3{-1.0 / 2, 3.0 / 2, 0}, zc),
        pair(L * Vec3{-1.0 / 2, -3.0 / 2, 0}, r * zc),
    };
    p.vertexPolyIsTerminals = {false, true, true, true};
    fillStraightPaths(p);
    p.validate();
    return p;
}

// Trifurcation: 2 junctions (vertices 1,2) + 4 terminals (vertices 3,4,5,6).
// Analogue of the Mathematica "5C_trifur" preset.  The two junctions are 3-way
// (each meets three edges) exactly like the bifurcation's single junction, so
// this mirrors that proven structure: a central junction-to-junction edge plus
// two spread-out terminals branching off each junction (a symmetric double-Y).
// Cross-sections are vertical bigons (z offset) and twists are zero, matching
// the bifurcation that converges to round tubes.
Preset makeTrifurcation(double lengthScale) {
    Preset p;
    p.name = "trifurcation";
    p.nVertices = 6;
    p.edgeList = {{1, 2}, {1, 3}, {1, 4}, {2, 5}, {2, 6}};
    p.edgeWidthList = {2 * PI, 2 * PI, 2 * PI, 2 * PI, 2 * PI};
    p.edgeTwistList = {0.0, 0.0, 0.0, 0.0, 0.0};

    p.vertexPolyLines = {{{1, 2}, {1, 2}, {1, 2}}, {{1, 2}, {1, 2}, {1, 2}},
                         {{1, 2}, {1, 2}},         {{1, 2}, {1, 2}},
                         {{1, 2}, {1, 2}},         {{1, 2}, {1, 2}}};
    p.vertexPolyLineLengths = {{PI, PI, PI}, {PI, PI, PI}, {PI, PI}, {PI, PI}, {PI, PI}, {PI, PI}};
    p.vertexPolyFacesByLines = {{{1, -2}, {2, -3}, {3, -1}}, {{1, -2}, {2, -3}, {3, -1}},
                                {{2, -1}},                   {{2, -1}},
                                {{2, -1}},                   {{2, -1}}};
    p.vertexPolyConnectionsByEdgeN = {{1, 2, 3}, {-1, 4, 5}, {-2}, {-3}, {-4}, {-5}};

    // In-plane layout (z = 0), scaled by lengthScale; vertical bigon cross-section.
    const double s = lengthScale;
    const Vec3 zc{0, 0, PI / 2};
    auto pair = [](const Vec3& base, const Vec3& off) {
        return std::vector<Vec3>{base + off, base - off};
    };
    p.vertexPolyCoordinates = {
        pair(s * Vec3{-1.5, 0.0, 0.0}, zc),  // junction 1
        pair(s * Vec3{1.5, 0.0, 0.0}, zc),   // junction 2
        pair(s * Vec3{-4.0, 2.5, 0.0}, zc),  // terminal off junction 1
        pair(s * Vec3{-4.0, -2.5, 0.0}, zc), // terminal off junction 1
        pair(s * Vec3{4.0, 2.5, 0.0}, zc),   // terminal off junction 2
        pair(s * Vec3{4.0, -2.5, 0.0}, zc),  // terminal off junction 2
    };
    p.vertexPolyIsTerminals = {false, false, true, true, true, true};

    fillStraightPaths(p);
    // Edge "lengths" follow the straight-line distance between the connected
    // vertex centers (sets each tube's mesh resolution along its axis).
    p.edgeLengthList.resize(p.numEdges());
    for (int e = 0; e < p.numEdges(); ++e) {
        const Vec3 a = meanCoords(p.vertexPolyCoordinates[p.edgeList[e][0] - 1]);
        const Vec3 b = meanCoords(p.vertexPolyCoordinates[p.edgeList[e][1] - 1]);
        p.edgeLengthList[e] = norm(b - a);
    }
    p.validate();
    return p;
}

}  // namespace msn
