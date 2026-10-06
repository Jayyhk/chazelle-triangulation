# Chazelle's Linear-Time Polygon Triangulation

A C++20 implementation of Chazelle's 1991 deterministic
linear-time algorithm for triangulating simple polygons.

See [DEVIATIONS.md](docs/DEVIATIONS.md) for confirmed departures from the cited
papers.

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

On Ubuntu/Debian, install the GMP development package before configuring:

```bash
sudo apt install libgmp-dev
```

## Usage

```bash
echo "n
x0 y0
x1 y1
...
x_{n-1} y_{n-1}" | ./build/chazelle
```

**Input:** vertex count n on the first line, then n lines of `x y` coordinates
for a simple polygon in either clockwise or counterclockwise boundary order.
Omit the repeated closing vertex. Decimal coordinates, including scientific
notation, are read exactly. Equal heights, horizontal edges, and consecutive
collinear vertices are supported.

**Output:** n−2 on the first line, followed by one line per triangle giving
three zero-based indices into the input vertex list, in clockwise order.

## Visualize

Add `--visualize` to save a picture of the triangulation:

```bash
echo "n
x0 y0
x1 y1
...
x_{n-1} y_{n-1}" | ./build/chazelle --visualize
```

Images are saved in `media/images/triangulation-YYYY-MM-DD.svg`.
Use `--visualize example.svg` to choose a filename.

## Animate

Install Manim with Python 3.11 or newer:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r src/animation/requirements.txt
```

See [Manim's installation instructions](https://docs.manim.community/en/stable/installation.html)
for native prerequisites.

Add `--animate` to render the algorithm:

```bash
echo "n
x0 y0
x1 y1
...
x_{n-1} y_{n-1}" | ./build/chazelle --animate
```

Videos are saved in `media/animations/algorithm-YYYY-MM-DD.mp4`.
Use `--animate example.mp4` to choose a filename.

## References

- **[C91]** B. Chazelle, "Triangulating a Simple Polygon in Linear Time,"
  _Discrete & Computational Geometry_, 6(5):485-524, 1991.
- **[FM84]** A. Fournier and D.Y. Montuno, "Triangulating Simple Polygons and
  Equivalent Problems," _ACM Transactions on Graphics_, 3(2):153-174, 1984.
- **[LT79]** R.J. Lipton and R.E. Tarjan, "A Separator Theorem for Planar
  Graphs," _SIAM Journal on Applied Mathematics_, 36(2):177-189, 1979.
