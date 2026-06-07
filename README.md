# min-surf-netw-cpp

A dependency-free C++17 port of the [min-surf-netw](https://github.com/x-meng-group/min-surf-netw)
Mathematica framework, which constructs and optimizes the geometry of physical
**network manifolds**. A network is modeled as a continuous surface made of
quad-meshed tubes (edges) that meet at polygonal junctions (vertices). The solver
finds an equilibrium 3D embedding by minimizing a tunable energy functional.

This port reimplements the algorithm in `QuadMeshNetworkScript.wl` function for
function, replacing the Mathematica built-ins it relied on:

| Mathematica | C++ replacement |
| --- | --- |
| `FindMinimum` (quasi-Newton) | L-BFGS with backtracking line search (`src/lbfgs.hpp`) |
| symbolic energy + analytic gradient | reverse-mode autodiff over a static tape (`src/autodiff.hpp`) |
| `ToElementMesh` cylinder FEM | direct structured quad grid (`Network::setupFEM`) |
| `Interpolation[..., Order->1]` | piecewise-linear interpolation (`src/interpolation.hpp`) |
| `RotationMatrix[{a,b}]` | analytic Rodrigues rotation (`src/linalg.hpp`) |
| `.mx` presets | JSON preset format + built-in generators (`src/preset.*`) |

There are **no external dependencies** — only the C++17 standard library.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This produces `build/min-surf-netw`.

## Run

```bash
# Built-in bifurcation preset (1 junction + 3 terminals)
./build/min-surf-netw --preset bifurcation --out bif

# Built-in trifurcation preset (2 junctions + 4 terminals)
./build/min-surf-netw --preset trifurcation --out tri

# From a JSON preset file
./build/min-surf-netw --json examples/bifurcation.json --out bif
```

Each run writes `<prefix>.obj` (one quad-mesh tube per edge, viewable in any
mesh viewer such as MeshLab or Blender) and `<prefix>_params.txt` (per-edge mesh
dimensions and the optimized `sqrtZScale`).

### Key options

| Option | Meaning | Default |
| --- | --- | --- |
| `--preset <name>` | built-in preset: `bifurcation` or `trifurcation` | `bifurcation` |
| `--json <file>` | load a preset from JSON instead of a built-in | — |
| `--scale <s>` | global FEM resolution (larger = finer mesh) | `1.0` |
| `--wfair <w>` | fairness weight (annealed down each stage) | `100` |
| `--wsurface <w>` | surface-area weight | `1e-2` |
| `--wcircle <w>` | terminal-circle weight | `1e3` |
| `--wpath <w>` | path-adherence weight | `0` |
| `--witer <n>` | number of `wfair` annealing stages (`/10` each) | `4` |
| `--maxiter <n>` | L-BFGS iterations per stage | `5000` |
| `--dump-preset` | write the chosen built-in preset to JSON and exit | — |

## Pipeline

The `Network` class mirrors the original function names and call order:

1. `setupFEM(scale)` — build a structured quad-cylinder mesh for every edge.
2. `doAttribute()` — path interpolation, parallel-transport frames, boundary rings.
3. `doConstraint()` — assemble the structural data for every energy term.
4. `doInitialize()` — place an initial 3D embedding along the guide paths.
5. `doFindMinimum(weights)` — minimize the energy with L-BFGS.

The energy functional combines (matching the Mathematica source): tube isometry,
interior and boundary fairness, soft boundary gluing between adjacent tubes,
terminal positioning constraints, and a surface-area term.

## Preset / JSON format

The original presets were Mathematica binary `.mx` files. Those are opaque, so
the C++ port defines an equivalent **JSON** schema that transcribes the same
data the `.mx` files would supply (`edgeList`, `vertexPolyLines`,
`vertexPolyFacesByLines`, `vertexPolyCoordinates`, …), keeping the original
1-based, signed-index conventions. The built-in `bifurcation` and `trifurcation`
generators emit this schema; see `examples/` for dumped instances. Because the
`.mx` contents are not available, the built-in presets are best-effort
reconstructions of the canonical topologies rather than byte-exact copies.

## License

GPL-3.0, matching the original project. See [LICENSE](LICENSE).
