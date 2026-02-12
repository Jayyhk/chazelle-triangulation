# Chazelle O(n) Simple Polygon Triangulation — Implementation Plan

## Overview

A C++20 implementation of Chazelle's 1991 linear-time simple polygon triangulation.
The algorithm computes a **visibility map** (horizontal trapezoidal decomposition) in
O(n) via a two-phase (up/down) submap refinement scheme, then converts it to a
triangulation in O(n) using Fournier & Montuno's Algorithms 2+3. The ray-shooting
oracle internally uses Lipton & Tarjan's O(n) planar separator theorem.

Every component is strictly O(n). No O(n log n) sorting or balanced-tree operations
are applied to n-sized inputs. All super-constant per-element costs occur on
submap-sized objects (O(n/γ) regions) whose totals telescope via geometric series.

### Dependency Graph

```
Lipton-Tarjan 1979              Fournier-Montuno 1984
(planar separator, O(n))        (trapezoid → triangulation, O(n))
         \                               /
          \                             /
           v                           v
             Chazelle 1991
      (visibility map in O(n))
```

### Key Parameters

| Symbol | Value       | Meaning                                                     |
| ------ | ----------- | ----------------------------------------------------------- |
| β      | 1/5         | Granularity exponent; canonical submaps are 2^⌈βλ⌉-granular |
| p      | ⌈log₂(n−1)⌉ | Number of grades; polygon padded to n = 2^p + 1             |
| γ      | 2^⌈βλ⌉      | Granularity at grade λ                                      |

### Oracle Cost Functions

These functions appear throughout the merge complexity bounds (Lemmas 3.1–3.5).
With β = 1/5 they evaluate to sub-linear powers of γ, ensuring the per-grade cost
decks geometrically.

| Function | Definition                                                     | Value (β = 1/5)                  | Role                                              |
| -------- | -------------------------------------------------------------- | -------------------------------- | ------------------------------------------------- |
| f(γ)     | Ray-shooting query cost in a γ-granular submap                 | O(γ^(11/15) · log γ) ≈ O(x^0.74) | Dominant oracle cost in fusion & conformality     |
| g(γ)     | Arc-cutting oracle: number of chain pieces per arc             | O(log γ) = O(βλ)                 | Multiplicative factor in conformality restoration |
| h(γ)     | Sub-granularity: granularity of sub-chains used by arc-cutting | O(γ^β) = O(γ^0.2)                | Additive term in conformality per-step cost       |

The key constraint from §4.1 is that the per-grade cost must satisfy:

> f(γ) · g(γ) · (h(γ) + log γ) = O(γ^(1−ε)) for some ε > 0

With β = 1/5 this yields a per-grade cost of O(n · 2^(−λ/76)), which sums
geometrically to O(n) across all p+1 grades.

---

## Project Layout

```
chazelle-triangulation/
├── CMakeLists.txt
├── PLAN.md
├── papers/
│   ├── chazelle1991-transcribed.tex
│   ├── fournier-montuno1984-transcribed.tex
│   └── lipton-tarjan1979-transcribed.tex
├── src/
│   ├── geometry/        # Step 1 — primitives
│   │   ├── point.h
│   │   ├── edge.h
│   │   ├── polygon.h
│   │   ├── polygon.cpp
│   │   └── perturbation.h
│   ├── separator/       # Step 2 — Lipton-Tarjan
│   │   ├── planar_graph.h
│   │   ├── planar_graph.cpp
│   │   ├── separator.h
│   │   └── separator.cpp
│   ├── visibility/      # Step 3 — submap data structures
│   │   ├── arc_structure.h
│   │   ├── chord.h
│   │   ├── submap.h         # Unified: tree, chords, arcs, tree decomp
│   │   ├── submap.cpp
│   │   ├── tree_decomposition.h
│   │   └── tree_decomposition.cpp
│   ├── oracles/         # Step 4 — ray-shooting & arc-cutting
│   │   ├── ray_shooting.h
│   │   ├── ray_shooting.cpp
│   │   ├── arc_cutting.h    # Thin utility (~10 lines of logic)
│   │   └── arc_cutting.cpp
│   ├── chazelle/        # Steps 5–6 — up-phase & down-phase
│   │   ├── fusion.h
│   │   ├── fusion.cpp
│   │   ├── conformality.h   # Callable from both up-phase and down-phase
│   │   ├── conformality.cpp
│   │   ├── granularity.h
│   │   ├── granularity.cpp
│   │   ├── merge.h
│   │   ├── merge.cpp
│   │   ├── grade_storage.h  # Per-grade canonical submaps + oracle structures
│   │   ├── up_phase.h
│   │   ├── up_phase.cpp
│   │   ├── down_phase.h
│   │   └── down_phase.cpp
│   ├── triangulate/     # Step 7 — Fournier-Montuno
│   │   ├── trapezoid.h      # Trapezoid type (only used in this step)
│   │   ├── vertex_node.h    # Doubly-linked polygon vertex list
│   │   ├── monotone.h
│   │   ├── monotone.cpp
│   │   ├── convert.h        # V(P) submap → trapezoid list conversion
│   │   ├── convert.cpp
│   │   ├── triangulate.h
│   │   └── triangulate.cpp
│   └── main.cpp         # Entry point: polygon in → triangles out
└── tests/               # Deferred
    ├── test_separator.cpp
    ├── test_visibility.cpp
    ├── test_oracles.cpp
    ├── test_merge.cpp
    ├── test_triangulate.cpp
    └── test_end_to_end.cpp
```

---

## Steps

### Step 1 — Core Geometry & Polygon Primitives

**Files:** `src/geometry/*`

**What to build:**

- `Point`: 2D point with double (or exact-arithmetic) coordinates.
- `Edge`: directed segment between two `Point`s, with index into the input table.
- `Polygon`: vertex list in boundary order, stored as the **read-only input table**
  (array of edges in boundary order, per Chazelle §2). If the polygon is closed,
  conceptually puncture a small hole to make it an open curve.
- **Padding:** pad vertex count to n = 2^p + 1 by inserting collinear vertices on
  the last edge. This is required so grades decompose cleanly into dyadic chains.
- **Symbolic y-perturbation** (`perturbation.h`): ensure no two vertices share the
  same y-coordinate. Use lexicographic comparison (y, then x, then vertex index)
  rather than mutating coordinates. All geometric predicates must route through this.
- **Double boundary (conceptual — NOT an explicit data structure):** The paper
  explicitly states (§2.3): "the notion of double boundary need not be encoded
  explicitly, i.e., no edges are duplicated in the table." The thickening of curve C
  into ∂C is handled entirely through **side flags** (`enum class Side { LEFT,
RIGHT }`) in arc-structures and geometric predicates. No `DoubleBoundary` class is
  constructed. The key properties:
  - **Regular vertices** (not local y-extrema, not endpoints): **2** companion
    vertices, one on each side of C.
  - **Local y-extrema** (not endpoints): **4** vertices — two companion pairs,
    one pair on each side. The "inside" pair produces a **null-length chord** and
    an **empty region**. These null-length chords and empty regions must be handled
    throughout the algorithm.
  - **Endpoints of C**: **2** companion/duplicate vertices, adjacent along ∂C and
    on both sides of C.
- **Double-backing arcs:** an arc can wrap around an endpoint of C, spanning both
  sides of the double boundary. This requires special-case logic in arc traversal
  and the arc-sequence table (§2.3).
- **Double identification:** finding which arc-structure corresponds to a given
  point on ∂C requires an O(log m) binary search in the arc-sequence table (§2.3).
  This operation is used throughout fusion, conformality restoration, and
  ray-shooting.
- **Spherical plane:** the algorithm operates on the spherical plane
  ([-∞,+∞]² with boundary identifications ≅ S²), not the Euclidean plane. Chords
  can wrap around infinity. The implementation must handle the unbounded exterior
  region correctly.

**No dependencies on other modules.**

---

### Step 2 — Lipton-Tarjan Planar Separator

**Files:** `src/separator/*`
**Paper reference:** Lipton & Tarjan 1979, §3 (Steps 1–10)

**What to build:**

- `PlanarGraph`: DCEL-style planar embedding. Each edge stores its two endpoints
  and four pointers (CW/CCW around each endpoint). Each vertex stores one incident
  edge. This is the input format required by the algorithm.
- `find_separator(graph, costs) → (A, B, C)`: the 10-step O(n) algorithm:
  1. Accept or compute planar embedding.
  2. Find connected components (Hopcroft-Tarjan DFS).
  3. BFS spanning tree; compute levels L(0), L(1), ….
  4. Find level l₁ such that cost of levels 0…l₁ ≤ ½.
  5. Find levels l₀ ≤ l₁ and l₂ ≥ l₁ with small total level size.
  6. Shrink levels l₀…l₂ into a single vertex; produce reduced graph.
  7. BFS tree of reduced graph; **triangulate all faces** of the new graph by
     adding non-tree edges as necessary. This face triangulation is a non-trivial
     sub-step: for each face, walk the facial cycle and add diagonals to make it
     a triangle. Required for step 9's cycle-replacement to work correctly.
  8. Find initial non-tree edge; compute fundamental cycle costs.
  9. Iterate: replace cycle edge to improve balance (each step removes a face;
     amortized O(n)).
  10. Construct final partition (A, B, C).
- `iterated_separator(graph, costs) → partition`: recursively apply
  `find_separator` to produce pieces of size ≤ μ^(2/3). Chazelle uses this in the
  ray-shooting oracle (§3.4) on the dual graph of S*, yielding separator set D* of
  size O(μ^(2/3)) and sub-pieces D₁, D₂, ….

**Depends on:** Step 1 (geometry types, only for vertex coordinates if embedding
needs to be computed).

---

### Step 3 — Visibility Submap Data Structures

**Files:** `src/visibility/*`
**Paper reference:** Chazelle 1991, §2

**What to build:**

- `Chord`: a horizontal segment connecting two mutually visible points on ∂C
  (the double boundary). A chord of V(C) partitions the double boundary cleanly
  (Lemma 2.1). A chord record stores (§2.3):
  - Pointers to the **2, 3, or 4 adjacent arc-structures** (not just 2 — there
    can be 3 or 4 arcs adjacent to a chord when chord endpoints coincide with
    vertices that produce extra arcs).
  - The chord IS a tree edge of the submap tree; the two regions it separates are
    the tree nodes at its endpoints — no separate region pointers needed.
  - No explicit y-coordinate: the chord's position is implicit from the arc-
    structures' input-table references.
- `ArcStructure`: encodes an arc (a connected piece of ∂C between two consecutive
  chord endpoints within a single region). Per §2.3:
  - **Single-edge arc** (t = 1): ONE pointer into the input table (the edge of P
    containing the arc) + ONE side flag (`enum class Side { LEFT, RIGHT }`).
  - **Multi-edge arc** (t > 1): TWO pointers into the input table (the edges
    containing the first and last sub-edges of the arc, in clockwise order) + TWO
    side flags.
  - A **back-pointer to the tree node** of the region this arc belongs to.
  - Null-length arcs must be representable (they occur at local extrema).
  - Endpoints of the arc are NOT stored — chords take care of that (§2.3).
- `Submap`: the central type. A submap is a coarsened visibility map obtained by
  retaining a subset of chords. **"Normal form" is a representation mode of this
  type, not a separate class.** A `Submap` in normal form contains (§2.3):
  1. **Submap tree** (unrooted): nodes = regions, edges = chords. Stored in
     **standard edge/node adjacency fashion** (adjacency lists, NOT parent/child
     pointers). The paper is explicit: "standard edge/node adjacency fashion."
  2. **Chord records**: each tree edge stores a chord record with pointers to 2–4
     adjacent arc-structures (see above).
  3. **Arc-sequence table**: arc-structures listed in canonical double-boundary
     traversal order.
  4. **Tree decomposition** (if conformal): stored separately as a rooted binary
     tree (see `TreeDecomposition` below).
     Key properties:
  - **Conformal**: submap tree has degree ≤ 4.
  - **γ-semigranular**: every region has weight **at most γ**.
  - **γ-granular**: γ-semigranular AND maximal — contracting any edge incident
    upon a node of degree < 3 would produce a node with weight **exceeding γ**.
    This is a maximality condition, not a per-region lower bound.
  - **Canonical**: conformal + γ-granular + normal-form.
  - **Weight** of a region: 0 if empty, otherwise the **maximum** number of
    non-null-length edges **in** (composing) any single arc of the region.
    Computing weight requires inspecting the arc-structures and comparing their
    input-table pointer spans (O(1) per arc).
  - The complete visibility map V(C) is itself a `Submap` — the maximal one with
    all chords retained. Use `using VisibilityMap = Submap;` as a type alias if
    desired, but no separate class is needed.
  - Region enumeration: provide an iterator interface over regions (tree nodes)
    that runs in O(r) where r = number of regions, not O(n).
- `TreeDecomposition`: hierarchical centroid decomposition of a conformal submap
  tree. **This IS rooted** (centroid-based, with parent/child pointers), unlike
  the submap tree itself which is unrooted. Internal nodes correspond to exit
  chords (centroid edges); leaves correspond to regions. Depth O(log r) where
  r is the number of regions (= number of chords + 1; for a γ-granular submap
  of an m-vertex curve, r = O(m/γ + 1) by Lemma 2.3). Computed in O(r log r)
  time (acceptable since r = O(m/γ) and the totals telescope). Used by
  conformality restoration (§3.2) for binary search over tree levels.

**Depends on:** Step 1.

---

### Step 4 — Oracles: Ray-Shooting & Arc-Cutting

**Files:** `src/oracles/*`
**Paper reference:** Chazelle 1991, §3.4

#### Ray-Shooting Oracle (Lemma 3.6)

**Purpose:** Given a horizontal ray from a point on ∂C, find the first chord or arc
of a submap S that the ray hits.

**Construction (preprocessing):**

1. Build **S\***, which is NOT the same object as S. S\* is a derived planar
   subdivision obtained by collapsing the double boundary to zero thickness and
   making all vertices (of C and chord endpoints) explicit on both sides. As a
   result, edges of S\* may be shorter than those of S, but unlike S, no edge of
   S\* has zero length (zero-length edges are contracted into vertices). Each face
   of S\* corresponds to exactly one region of S, though faces may have more
   vertices incident upon them.
2. Compute the **dual graph G** of S\*. G is planar with μ = O(m/γ + 1) nodes
   (one per face of S\*). Two faces of S\* are adjacent if and only if either
   (a) they share a chord, or (b) one of them has a chord endpoint that abuts on
   a non-null-length arc of the other's region. Adjacencies of type (a) come
   directly from the submap tree edges. Adjacencies of type (b) require double
   identification (§2.3) — locating which arc passes through a chord endpoint.
   **G is used as the input to the Lipton-Tarjan separator, NOT S\* directly.**
3. Apply iterated Lipton-Tarjan separator (Step 2) to G, producing:
   - Separator set D\* of size O(μ^(2/3)) (a set of nodes/faces).
   - Sub-pieces D₁, D₂, … each of size ≤ μ^(2/3), with region-to-piece
     membership lookups. This is a **hierarchical partition** of the dual graph
     nodes, not a flat partition.
4. **Vertical line structure:** Take a vertical line passing to the right of all
   vertices of P, and intersect it with the chords of the regions in S. This
   breaks up the line into segments, each falling entirely within some region.
   Store this as a sorted list for O(log μ) binary search — used when a ray
   doesn't hit any separator region and needs initial region identification.

**Query:** O(γ^(1/3) · m^(2/3)) time.

- Walk the separator hierarchy; at each level, test the separator faces (local
  shooting = check ≤ 4 arcs per region, since conformality guarantees degree ≤ 4).
- Descend into the appropriate sub-piece; repeat.

**Preprocessing cost:** O(m log m / γ + 1) per submap. Telescopes to O(n) across
all grades.

#### Arc-Cutting Oracle

**Purpose:** Given a subarc of ∂C, decompose it into O(log γ) chains from prior
grades whose canonical submaps have already been computed (in the up-phase).

**Implementation:** This is a **thin utility function** (~10 lines of logic), not a
heavyweight module. It performs binary decomposition of the arc's vertex range into
dyadic intervals, each matching a chain from some grade λ' < λ. Returns pointers to
the stored canonical submaps for these chains from `grade_storage`.

Concretely: take the arc's endpoint positions in the input table, decompose the
index range into O(log γ) dyadic intervals, look up each interval's corresponding
chain and its precomputed canonical submap from `grade_storage[λ'][chain_idx]`.

**Cost:** O(log γ) = O(βλ) pieces per arc. Each piece has a precomputed canonical
submap and ray-shooting structure.

**Depends on:** Steps 1, 2, 3.

---

### Step 5 — Up-Phase: Three-Stage Merge & Grade Processing

**Files:** `src/chazelle/fusion.*`, `src/chazelle/conformality.*`,
`src/chazelle/granularity.*`, `src/chazelle/merge.*`, `src/chazelle/up_phase.*`
**Paper reference:** Chazelle 1991, §3 (merge) and §4.1 (up-phase)

#### Three-Stage Merge (§3)

Given two canonical submaps S₁ (of chain C₁, n₁ vertices, γ₁-granular) and S₂
(of chain C₂, n₂ vertices, γ₂-granular) where C = C₁ ∪ C₂:

**Stage 1 — Fusion (§3.1):**

- Fusion consists of **two symmetric passes**: first fuse S₁ into S₂ (walk S₁'s
  exit-chord endpoints while querying S₂), then fuse S₂ into S₁ (walk S₂'s
  exit-chord endpoints while querying S₁). These are NOT a single interleaved walk.
- **Pass 1 (fuse S₁ into S₂):** Walk exit-chord endpoints a₀, a₁, …, a\_{m+1} of
  S₁ along ∂C₁. For each aᵢ, determine what it sees horizontally by querying S₂
  via ray-shooting. This discovers new chords crossing between the two chains.
- **Start-up:** find what a₀ sees (initial ray-shoot into S₂).
- **Main loop:** advance pointer p along ∂C₂, maintaining invariant that p tracks
  the arc of S₂ currently visible. At each step, either advance to the next S₁
  chord endpoint or discover a new chord.
- **Junction vertex C₁ ∩ C₂:** the companion vertices a₀ and a\_{m+1} from the
  duplication of the shared vertex require special-case handling, especially when
  the junction is a local extremum (producing 4 double-boundary vertices).
- **Pass 2 (fuse S₂ into S₁):** repeat symmetrically.
- Sort the O(n₁/γ₁ + n₂/γ₂) discovered chord endpoints along ∂C (by input-table
  edge index, then y-coordinate). This sort is on a submap-sized set, not n.
- Output: the **fusion submap** S — a valid submap of V(C) but not yet conformal or
  appropriately granular.
- **Cost:** O((n₁/γ₁ + n₂/γ₂ + 1)(f(γ₂) + log(n₁+n₂))) per Lemma 3.1.

**Stage 2 — Conformality Restoration (§3.2):**

- The fusion submap S may have regions with > 4 arcs (tree nodes with degree > 4).
- For each such region R, Lemma 3.3 guarantees two non-consecutive arcs Aᵢ, Aⱼ
  with a vertex of ∂C on Aᵢ visible to Aⱼ.
- Find such a visible pair via the **arc-cutting oracle**: decompose the
  problematic arc into O(g(γ)) subarcs, each with a canonical submap Sα from a
  prior grade. Then **binary search through the tree decomposition of Sα**: at each
  level, test O(1) diametrical chords using ray-shooting and the topological lemma
  (Lemma 2.4) to narrow down the correct sub-region.
- Insert the discovered chord, splitting R into two regions each with fewer arcs.
- Repeat until every region has ≤ 4 arcs (i.e., tree degree ≤ 4).
- **Cost:** O((n₁/γ₁ + n₂/γ₂ + 1) · f·g·(h + log γ₂)) per Lemma 3.4.

**Stage 3 — Granularity Enforcement (§3.3):**

- The conformal submap is already **γ-semigranular** (every region has weight
  ≤ γ) because conformality restoration (§3.2) did not remove any exit chords,
  and no arc has more than γ₂ edges. What may fail is the **maximality
  condition** (condition (ii) of γ-granularity).
- §3.3 **only removes chords** — it never inserts any. For each exit chord,
  check whether contracting the corresponding tree edge (merging the two
  adjacent regions) would produce a node with weight ≤ γ. If so, the chord is
  redundant: remove it. Chords need be processed only once, since removals
  cannot make a previously non-removable chord removable.
- This is a single pass over the submap tree.
- **Cost:** O(n₁/γ₁ + n₂/γ₂ + 1) — just tree traversal.

**Result:** a canonical (conformal, γ-granular, normal-form) submap of V(C).

#### Grade Processing (§4.1)

- Grades λ = 0, 1, …, p. At grade λ, there are (n−1)/2^λ chains, each with
  2^λ + 1 vertices.
- Grade 0: each chain is a single edge; its visibility map and canonical submap are
  trivial.
- Grade λ > 0: each chain C at grade λ is the union of two grade-(λ−1) chains C₁,
  C₂. Merge their canonical submaps (from grade λ−1) using the three-stage merge.
  Then preprocess the resulting canonical submap for ray-shooting (Step 4).
- **Grade storage** (`grade_storage.h`): all canonical submaps and ray-shooting
  structures must be retained for the down-phase and arc-cutting oracle.
  Organization: `std::vector<std::vector<CanonicalSubmap>> grade_storage(p+1)`
  where `grade_storage[λ]` has `(n−1)/2^λ` entries. Each entry stores the
  canonical submap plus its preprocessed ray-shooting structure. **Total storage
  is O(n)** by geometric series: at grade λ, total regions across all chains is
  O(n/2^(βλ) + n/2^λ), which sums to O(n) over all λ.
- **Cost per grade:** O(n · 2^(−λ/76)) (with β = 1/5).
- **Total across all grades:** geometric series → O(n).

**Depends on:** Steps 1, 2, 3, 4.

---

### Step 6 — Down-Phase: Recursive Refinement

**Files:** `src/chazelle/down_phase.*`
**Paper reference:** Chazelle 1991, §4.2

**Input:** the canonical submap of V(P) from the top grade (λ = p) of the up-phase,
plus the **entire `grade_storage`** array — all canonical submaps and ray-shooting
structures stored during the up-phase for every chain at every grade. The paper is
explicit (§4.2): "it uses data structures left behind during the [up-]phase (the
submaps for the chains and their ray-shooting structures)."

**Algorithm:**

- Process grades λ = p, p−1, …, 1.
- At each grade, refine the current submap by filling in details within each region:
  1. For each region R of the current submap, enumerate its arcs (≤ 4 by
     conformality).
  2. **Symbolic tilting:** the exit chords bounding R are temporarily treated as
     edges of the input curve (not as chords) by symbolically tilting them. After
     merging, they are reinterpreted as chords, which may break conformality.
     This requires a **full §3.2 conformality restoration pass** (the same
     procedure used during the up-phase merge, NOT a simplified version).
     `conformality.cpp` must be designed to be callable from both the up-phase
     and the down-phase.
  3. Use the **arc-cutting oracle** to decompose each arc into O(log γ) chains from
     prior grades.
  4. Retrieve the canonical submaps for these chains (stored during the up-phase).
  5. Merge these submaps into a **2^⌈β⌈βλ⌉⌉-granular conformal submap** of V(R\*)
     — NOT the full visibility map. This is a coarser approximation.
  6. Extract the exit chords falling within each region R; add them to refine the
     global submap.
  7. **Recurse** by reducing the granularity parameter from λ to ⌈βλ⌉ (induction
     on λ, while l — the grade of the chain — stays fixed). Since ⌈βλ⌉ ≤ λ−1
     for large λ, this terminates. **Base case:** λ is a small constant, so
     regions have bounded size and can be filled in directly in O(1) per region.
- At grade 0, every region contains at most one edge of P; the visibility map is
  complete.

**Complexity:** Lemma 4.2 shows the recurrence solves to (c − 1/λ) · 2^l per chain
of size 2^l at grade λ. Summed across all chains and grades, this telescopes to
O(n).

**Output:** the complete visibility map V(P) as a normal-form `Submap` with all
chords retained (every region is a single trapezoid, i.e., granularity 1). This is
converted to the Fournier-Montuno trapezoid list format by Step 7's bridging step.

**Depends on:** Steps 1, 3, 4, 5.

---

### Step 7 — Fournier-Montuno: Visibility Map → Triangulation

**Files:** `src/triangulate/*`
**Paper reference:** Fournier & Montuno 1984, Algorithms 2 & 3 (NOT Algorithm 1)

**Bridging step — V(P) submap → trapezoid list** (`convert.h/.cpp`):

The down-phase outputs V(P) as a normal-form `Submap` with all chords retained
(every region is a single trapezoid). Fournier-Montuno needs a list of `Trapezoid`
objects with explicit top/bottom vertices and left/right edge references, plus a
doubly-linked vertex list. This conversion is a **linear-time traversal** of the
submap tree:

- For each region (tree node), extract the bounding chord y-coordinates (top/bottom
  vertices) and the left/right polygon edges from the arc-structures.
- Build a `Trapezoid` for each region: top vertex, bottom vertex, left edge ref,
  right edge ref.
- Build the doubly-linked `VertexNode` list in polygon boundary order.

**Data structures for Algorithms 2+3:**

- `Trapezoid`: top vertex, bottom vertex, left edge reference, right edge reference.
  Only used in this step (not in Chazelle's submap machinery).
- `VertexNode`: intrusive doubly-linked list node for polygon vertices. Contains
  prev/next pointers, a pointer to the associated trapezoid (`trapezoid(vertex)`),
  and the vertex coordinates. Algorithm 2 rewires prev/next pointers to create
  sub-polygons; Algorithm 3 removes vertices from this list.

**Input:** the complete visibility map V(P) from Step 6, converted to the above
format.

**Algorithm 2 — Trapezoid classification & diagonal insertion:**

- Each trapezoid has two defining vertices (top and bottom).
- **Class A:** the two defining vertices share a polygon edge → already part of the
  triangulation boundary.
- **Class B:** the two defining vertices do NOT share a polygon edge → insert a
  diagonal between them. This diagonal becomes a triangulation edge.
- After all Class B diagonals are inserted, the polygon is decomposed into
  **unimonotone subpolygons** (one monotone chain + one single edge).
- Algorithm 2 is **recursive**: each Class B diagonal splits the polygon into two
  sub-polygons via **pointer surgery** on the doubly-linked vertex list:
  `save_next = next(current); save_prev = prev(bottom);
next(current) = bottom; prev(bottom) = current;`
  Each sub-polygon is then processed recursively. When no Class B trapezoid remains,
  the piece is unimonotone and handed to Algorithm 3.
- **Cost:** O(n) total (amortized over the recursion tree: each diagonal insertion
  is O(1), there are O(n) diagonals, and Algorithm 3 costs O(n) total).

**Algorithm 3 — Unimonotone polygon triangulation:**

- **Start vertex:** one of the two endpoints of the base edge, chosen as the
  **topmost** if the monotone chain is on the right, or the **bottommost** if it
  is on the left. Set `current = next(start)` to begin walking the monotone chain.
- Walk the monotone chain. At each vertex, check if
  `angle(prev(current), current, next(current))` is **convex**.
  - If convex: remove the vertex (emit a triangle formed by prev, current, next),
    then backtrack one vertex (the predecessor may have become convex).
  - If reflex: advance to the next vertex.
- The backtracking is amortized O(1): total forward steps + backward steps ≤ 2n.
- **Cost:** Θ(n) per subpolygon; Θ(n) total across all subpolygons.

**Output:** a list of triangles (triples of vertex indices) forming a triangulation
of the input polygon. Total cost: O(n).

**Depends on:** Steps 1, 6.

---

## Build & Entry Point

**`CMakeLists.txt`:** C++20, single executable target linking all `src/` modules.
No external dependencies beyond the C++20 standard library.

**`src/main.cpp`:** reads a polygon (vertex list), runs the full pipeline
(pad → up-phase → down-phase → Fournier-Montuno), and outputs the triangle list.

---

## Implementation Order

| Phase | Modules                             | Rationale                                                                                           |
| ----- | ----------------------------------- | --------------------------------------------------------------------------------------------------- |
| 1     | Step 1 (geometry)                   | Foundation; no dependencies                                                                         |
| 2     | Step 2 (separator)                  | Independent algorithm; needed by oracles                                                            |
| 3     | Step 7 (Fournier-Montuno)           | Independent algorithm; final pipeline stage; testable independently on hand-built trapezoidizations |
| 4     | Step 3 (visibility data structures) | Central types for the Chazelle core                                                                 |
| 5     | Step 4 (oracles)                    | Requires separator + visibility types                                                               |
| 6     | Step 5 (up-phase)                   | Core algorithm, first half                                                                          |
| 7     | Step 6 (down-phase)                 | Core algorithm, second half                                                                         |
| 8     | Integration & main                  | Wire everything together                                                                            |

---

## Implementation-Critical Subtleties

These details from the papers are easy to miss but essential for correctness:

1. **Weight is a max, not a sum (§2.2):** A region's weight is 0 if empty, otherwise
   the **maximum** number of non-null-length polygon edges composing any single arc
   of that region. This is NOT the total edge count or the number of V(C)
   sub-regions.

2. **Null-length chords and empty regions (§2.1):** Local y-extrema produce
   null-length chords and empty (weight-0) regions. These propagate through every
   submap operation — fusion, conformality, granularity, and the down-phase must all
   handle them correctly.

3. **Double-backing arcs (§2.3):** An arc can wrap around an endpoint of C, spanning
   both sides of the double boundary. Arc traversal, the arc-sequence table, and
   double identification must account for this.

4. **Double identification (§2.3):** Locating which arc passes through a given point
   requires O(log m) binary search in the arc-sequence table. Used pervasively in
   fusion and conformality restoration.

5. **Spherical plane topology (§2):** Everything is embedded in [-∞,+∞]² ≅ S².
   Chords can extend to/from ±∞. The unbounded exterior is a single region.

6. **Junction vertex handling in fusion (§3.1):** The companion vertices a₀ and
   a\_{m+1} at the shared vertex C₁ ∩ C₂ need special-case logic, particularly when
   the junction is a local extremum (4 double-boundary vertices, null-length chord).

7. **Symbolic chord tilting in the down-phase (§4.2):** When processing region R,
   its bounding chords are temporarily reclassified as curve edges (symbolically
   tilted). After merging, they revert to chords, potentially breaking conformality.
   This requires the **full §3.2 conformality restoration procedure** — the same
   mechanism used during the up-phase merge (Lemma 3.2 binary search through tree
   decomposition, Lemma 3.3 existence guarantee). Design `conformality.cpp` to be
   callable from both the up-phase merge and the down-phase refinement.

8. **Lipton-Tarjan Step 9 — non-tree-edge sub-case (§3):** When neither (vᵢ, y) nor
   (y, wᵢ) is a tree edge, the algorithm must find the tree path from y to the
   current cycle, then scan edges inside both sub-cycles alternately. This is the
   most complex case in the separator algorithm and requires careful implementation.

---

## Testing (Deferred)

Unit tests are planned but not a priority. When implemented:

| Test file              | What it covers                                                                                                 |
| ---------------------- | -------------------------------------------------------------------------------------------------------------- |
| `test_separator.cpp`   | Lipton-Tarjan on known planar graphs (grids, trees, random); verify partition sizes and separator bound        |
| `test_visibility.cpp`  | Submap construction, conformality checks, normal-form invariants on small polygons                             |
| `test_oracles.cpp`     | Ray-shooting correctness against brute-force; arc-cutting decomposition validity                               |
| `test_merge.cpp`       | Three-stage merge on hand-crafted submap pairs; verify conformality + granularity post-merge                   |
| `test_triangulate.cpp` | Fournier-Montuno on hand-built trapezoidizations; verify triangle count = n−2                                  |
| `test_end_to_end.cpp`  | Full pipeline on convex polygons, monotone polygons, stars, random simple polygons; verify valid triangulation |

---

## Constants & Configuration

| Constant            | Value       | Description                                                                   |
| ------------------- | ----------- | ----------------------------------------------------------------------------- |
| `BETA`              | 1.0/5.0     | Granularity exponent (constexpr double); controls how fast γ grows with grade |
| `MAX_TREE_DEGREE`   | 4           | Conformality bound: max degree of any node in a submap tree                   |
| `granularity(λ)`    | 2^⌈βλ⌉      | Computed per grade; canonical submaps at grade λ are this-granular            |
| `GRANULARITY_UPPER` | γ           | Max weight per region (semigranular/granular condition (i): weight ≤ γ)       |
| `NUM_GRADES(n)`     | ⌈log₂(n−1)⌉ | Total number of grades; polygon padded to n = 2^p + 1                         |
