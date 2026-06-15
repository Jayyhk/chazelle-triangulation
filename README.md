# Chazelle's Linear-Time Polygon Triangulation

A faithful C++20 implementation of Chazelle's 1991 deterministic linear-time
algorithm for triangulating simple polygons [C91]. The implementation includes
the Lipton-Tarjan [LT79] planar separator for ray-shooting and the
Fournier-Montuno [FM84] method for the final trapezoidization-to-triangulation
step.

## Implementation Progress

**[C91]** Chazelle's Visibility Algorithm
- [x] §2 — Visibility Maps and Submaps
  - [x] §2.0 — Spherical Plane and Orientation
  - [x] §2.1 — The Visibility Map
  - [x] §2.2 — Visibility Submaps
  - [x] §2.3 — Conformality and Granularity
  - [x] §2.4 — Representation Issues
  - [x] §2.5 — A Topological Lemma
- [ ] §3 — Merging Two Submaps
  - [x] §3.0 — Merge Setup and Oracle Primitives
  - [ ] §3.1 — Fusion of Two Submaps (`merge.cpp::fuse` is still TODO)
  - [ ] §3.2 — Restoring Conformality
  - [ ] §3.3 — Maintaining Granularity
  - [ ] §3.4 — Implementing the Oracles
- [ ] §4 — The Visibility Algorithm
  - [ ] §4.0 — Chains, Grades, and Algorithm Overview
  - [ ] §4.1 — The Up-Phase
  - [ ] §4.2 — The Down-Phase

**[LT79]** Lipton-Tarjan Planar Separator
- [ ] §3 — An Algorithm for Finding a Good Partition

**[FM84]** Fournier-Montuno Triangulation
- [ ] Extract trapezoids from V(P)
- [ ] Algorithm 2 — Triangulating the Trapezoidized Polygon
- [ ] Algorithm 3 — Triangulating Unimonotone Polygons

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

## References

- **[C91]** B. Chazelle, "Triangulating a Simple Polygon in Linear Time,"
  _Discrete & Computational Geometry_, 6(5):485-524, 1991.
- **[FM84]** A. Fournier and D.Y. Montuno, "Triangulating Simple Polygons and
  Equivalent Problems," _ACM Transactions on Graphics_, 3(2):153-174, 1984.
- **[LT79]** R.J. Lipton and R.E. Tarjan, "A Separator Theorem for Planar
  Graphs," _SIAM Journal on Applied Mathematics_, 36(2):177-189, 1979.
