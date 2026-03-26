# Chazelle's Linear-Time Polygon Triangulation

A faithful C++20 implementation of Chazelle's 1991 deterministic linear-time
algorithm for triangulating simple polygons [C91]. The implementation includes
the Lipton-Tarjan [LT79] planar separator for ray-shooting and the
Fournier-Montuno [FM84] method for the final trapezoidization-to-triangulation
step. Intended as a reference implementation for studying the algorithm.

> [!CAUTION]
> This is a work in progress! With high probability there are bugs present. Use at your own risk.

## Algorithm

The pipeline has five stages:

1. **Pad** the polygon to 2^p + 1 vertices.
2. **Up-phase:** build canonical submaps bottom-up across grades [C91 &sect;4.1].
3. **Down-phase:** refine to the full visibility map V(P) top-down [C91 &sect;4.2].
4. **Convert** V(P) to a trapezoidization.
5. **Triangulate** monotone pieces [FM84].

Total: O(n) time.

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Requires a C++20 compiler (GCC 12+, Clang 15+) and CMake 3.20+.

## Usage

```bash
echo "N
x0 y0
x1 y1
...
x_{N-1} y_{N-1}" | ./build/chazelle
```

**Input:** vertex count N on the first line, then N lines of `x y` coordinates
in counterclockwise boundary order.

**Output:** timing information per phase, followed by triangle count and one
line per triangle giving three vertex indices.

## Project Structure

```
src/
├── polygon/         Polygon, edges, SoS perturbation              [C91 §2]
├── submap/          Submaps, arcs, chords, tree decomposition     [C91 §2.1–2.5]
├── merge/           Core merge operation + oracles                [C91 §3]
│   └── separator/   Lipton-Tarjan planar separator                [LT79]
├── phases/          Algorithm driver: up/down phases              [C91 §4]
├── triangulate/     Fournier-Montuno triangulation                [FM84]
├── common.h
└── main.cpp         Entry point
tests/               Unit and integration tests
papers/              Transcribed reference papers
visualizer/          Python visualization tools
```

## References

- **[C91]** B. Chazelle, "Triangulating a Simple Polygon in Linear Time,"
  _Discrete & Computational Geometry_, 6(5):485-524, 1991.
- **[FM84]** A. Fournier and D.Y. Montuno, "Triangulating Simple Polygons and
  Equivalent Problems," _ACM Transactions on Graphics_, 3(2):153-174, 1984.
- **[LT79]** R.J. Lipton and R.E. Tarjan, "A Separator Theorem for Planar
  Graphs," _SIAM Journal on Applied Mathematics_, 36(2):177-189, 1979.

## License

<!-- TODO: add license -->
