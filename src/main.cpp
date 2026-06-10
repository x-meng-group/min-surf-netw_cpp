// Command-line driver for the network-manifold solver.
//
// Usage:
//   min-surf-netw [options]
//
// Options:
//   --preset <name>     built-in preset: bifurcation | trifurcation  (default bifurcation)
//   --json <file>       load a preset from a JSON file instead of a built-in
//   --scale <s>         global FEM resolution scale (default 1.0; larger = finer)
//   --wfair <w>         fairness weight (default 100)
//   --wsurface <w>      surface-area weight (default 1e-2)
//   --wcircle <w>       terminal-circle weight (default 1e3)
//   --wpath <w>         path-mean weight (default 0)
//   --witer <n>         number of wfair annealing stages (default 4; wfair /= 10 each)
//   --maxiter <n>       L-BFGS max iterations per stage (default 5000)
//   --out <prefix>      output file prefix (default "out")
//   --dump-preset       write the chosen preset to <prefix>_preset.json and exit
//   --quiet             suppress per-iteration logging
#include <cstdio>
#include <cstring>
#include <string>

#include "io.hpp"
#include "network.hpp"
#include "preset.hpp"

using namespace msn;

int main(int argc, char** argv) {
    std::string presetName = "bifurcation";
    std::string jsonFile;
    std::string outPrefix = "out";
    double scale = 1.0;
    double wfair = 100.0, wSurface = 1e-2, wCircle = 1e3, wPath = 0.0;
    int wIterStages = 4;
    int maxIter = 5000;
    bool dumpPreset = false;
    bool quiet = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* name) -> std::string {
            if (i + 1 >= argc) { std::fprintf(stderr, "Missing value for %s\n", name); std::exit(2); }
            return argv[++i];
        };
        if (a == "--preset") presetName = next("--preset");
        else if (a == "--json") jsonFile = next("--json");
        else if (a == "--scale") scale = std::stod(next("--scale"));
        else if (a == "--wfair") wfair = std::stod(next("--wfair"));
        else if (a == "--wsurface") wSurface = std::stod(next("--wsurface"));
        else if (a == "--wcircle") wCircle = std::stod(next("--wcircle"));
        else if (a == "--wpath") wPath = std::stod(next("--wpath"));
        else if (a == "--witer") wIterStages = std::stoi(next("--witer"));
        else if (a == "--maxiter") maxIter = std::stoi(next("--maxiter"));
        else if (a == "--out") outPrefix = next("--out");
        else if (a == "--dump-preset") dumpPreset = true;
        else if (a == "--quiet") quiet = true;
        else if (a == "--help" || a == "-h") {
            std::printf(
                "Usage: %s [--preset bifurcation|trifurcation] [--json file] [--scale s]\n"
                "          [--wfair w] [--wsurface w] [--wcircle w] [--wpath w]\n"
                "          [--witer n] [--maxiter n] [--out prefix] [--dump-preset] [--quiet]\n",
                argv[0]);
            return 0;
        } else {
            std::fprintf(stderr, "Unknown argument: %s\n", a.c_str());
            return 2;
        }
    }

    Preset preset;
    try {
        if (!jsonFile.empty()) {
            preset = loadPresetFile(jsonFile);
        } else if (presetName == "bifurcation") {
            preset = makeBifurcation();
        } else if (presetName == "trifurcation") {
            preset = makeTrifurcation();
        } else {
            std::fprintf(stderr, "Unknown preset '%s' (use bifurcation or trifurcation)\n",
                         presetName.c_str());
            return 2;
        }
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "Failed to load preset: %s\n", ex.what());
        return 1;
    }

    if (dumpPreset) {
        std::string path = outPrefix + "_preset.json";
        std::ofstream f(path);
        f << presetToJson(preset);
        std::printf("Wrote preset to %s\n", path.c_str());
        return 0;
    }

    std::printf("Preset: %s  (%d vertices, %d edges)\n", preset.name.c_str(), preset.nVertices,
                preset.numEdges());

    try {
        Network net(preset);
        net.setupFEM(scale);
        net.doAttribute();
        net.doConstraint();
        net.doInitialize(1.0, 1.0);

        double w = wfair;
        for (int stage = 0; stage < wIterStages; ++stage) {
            SolveWeights sw;
            sw.wfair = w;
            sw.wSurface = wSurface;
            sw.wCircle = wCircle;
            sw.wPathPerFair = wPath;
            sw.maxIterations = maxIter;
            sw.verbose = !quiet;
            std::printf("[stage %d/%d] wfair=%.6g\n", stage + 1, wIterStages, w);
            net.doFindMinimum(sw);
            w /= 10.0;
        }

        std::string objPath = outPrefix + ".obj";
        std::string dumpPath = outPrefix + "_params.txt";
        writeObj(net, objPath);
        writeParamDump(net, dumpPath);
        std::printf("Wrote mesh to %s and parameters to %s\n", objPath.c_str(), dumpPath.c_str());
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "Solver error: %s\n", ex.what());
        return 1;
    }
    return 0;
}
