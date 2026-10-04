# Confirmed deviations from the cited papers

Only **established, current departures** belong here. Each entry identifies the
paper instruction, a concrete contradictory implementation behavior, and its
complexity consequences. Ambiguous wording, hypotheses, corrected historical
exceptions, and behavior required by another explicit paper instruction are
excluded. Absence from this list is not a certificate of compliance.

Scope covers the implemented C++ visibility-map construction through C91 §4.2
and its LT79 separator, plus FM84 visibility-to-trapezoid extraction and
Algorithms 2–3, including triangle output and vertex adjacency. `FM84 tex N`
refers to [FM84's local transcription](../papers/fournier-montuno1984-transcribed.tex).
`C91 tex N` and
`LT79 tex N` refer to [C91's local transcription](../papers/chazelle1991-transcribed.tex)
and [LT79's local transcription](../papers/lipton-tarjan1979-transcribed.tex).
Publication-error findings below were also checked against the original PDFs.

## 1. C91 §3.3: removable chords require bounded rechecks

**Paper:** C91 tex 276 says a removal cannot make a previously nonremovable
chord removable, and consequently each chord needs processing only once.
This statement is also in the original publication, p. 512 (PDF page 28).

**Established contradiction:** criterion (ii), tex 121, only tests contractions
incident upon a node of degree less than three. An edge joining two degree-3
nodes is initially ineligible. Contracting a leaf attached to either node
lowers that node's degree to two, making the surviving edge eligible. When γ
exceeds the entire doubled-boundary weight, it is then removable.
`tests/s3.3-tests.cpp::test_degree_drop_recheck` constructs an actual conformal
submap, orders its hub chords before its leaf chords, and checks complete
contraction. A once-only variant fails the criterion-(ii) postcondition.
The publication's arbitrary-order monotonicity assertion is false.

**Implemented correction:** after contraction, enqueue the surviving node's
at-most-four incident chords. This enforces the original granularity definition
and Lemma 3.5; it changes the once-only processing instruction.

**Time:** with m initial chords, there are m initial work items, at most m
removals, and at most 4m added work items. Conformality bounds each contraction
and test by O(1). Time and space remain O(m), preserving §3.3 and overall O(n).

## 2. C91 §§3.0/4.1: individual piece granularity differs from uniform h

**Paper:** tex 170 requires each non-boundary cutter piece to have an
h(γ_i)-granular submap. Tex 343's construction instead reuses canonical
submaps whose individual granularities are *at most* h. Both statements
appear in the original publication, pp. 501 and 517 (PDF pages 17 and 33).

**Established contradiction:** full granularity is not monotone in γ, because
criterion (ii), tex 121, requires more contraction when γ increases. The
production result checked by `tests/s4.1-tests.cpp::test_up_phase_cutter`, at
λ=3, includes full grade-0, vertex-to-vertex pieces with γ_j=1 and two wrap
chords. They are 1-granular but **not 2-granular**, while the cutter advertises
h=2. These are not partial-edge pieces exempted by tex 170. The regression
checks actual granularity, failure of uniform-h granularity, and h-semigranularity.
Thus the §4.1 reuse construction does not literally establish the stronger
§3.0 contract.

**Implemented correction:** `ArcPiece` declares γ_j≤h and validates full
granularity at γ_j. The search uses h as a bound on region weight; it does not
claim full h-granularity or contract pieces during a query.

**Time:** a piece has at most γ_i edges (tex 341). Lemma 2.3 at γ_j≥1 bounds
its regions by O(γ_i), giving O(log γ_i) centroid-tree depth. Conformality and
weight ≤γ_j≤h give O(h) leaf scanning. These are precisely the bounds used
by Lemma 3.2, tex 248. The reused maps keep cutting O(g), and the f/g/h and
overall O(n) bounds survive the corrected contract.

## 3. C91 §4.1: resetting inputs also rebuilds normal form

**Paper:** tex 339 claims initial resetting takes O(M) over the input submap
sizes. The original publication makes the same claim on p. 517 (PDF page 33).
Normal form includes the centroid decomposition (tex 139), and merge inputs
must be in normal form (tex 166).

**Established implementation departure:** `compute_canonical_portion` resets
its copied inputs through `reset_chain_granularity`, which calls
`enforce_granularity` and then `normalize`. That normalization
sorts arcs using `std::sort`. The centroid decomposition now takes O(M), but
the comparison sort still gives the local reset/normalization stage an
O(M log M) bound, rather than tex 339's O(M) bound. The general normalization
bound is explicitly permitted at tex 116 and tex 276. This establishes an
implementation efficiency difference; it does **not** prove that a linear
normalization algorithm is impossible or that the paper's bound cannot be
achieved by another representation.

**Overall time:** let β=1/5 and n be the original input vertex count. Minimal
padding at tex 316 keeps the processed curve below 2n vertices. Lemma 2.3
(tex 126) bounds a stored grade-k submap by
M_k=O(2^(k−⌈βk⌉)+1) records, including its tree decomposition. Tex 339's
partition uses at most two maps per grade k<λ. Copying, contracting, and
normalizing these maps therefore costs

`A_λ = O(sum_{k<λ} M_k log(M_k+1)) = O(λ·2^((1−β)λ))`.

For the up-phase, tex 319 and tex 355 give O(n·2^(−λ)) chains at grade λ.
The aggregate reset cost is consequently

`O(n·sum_{λ>=1} λ·2^(−βλ)) = O(n)`.

For the down-phase, `RegionBoundaryGeometry::canonical_original_piece` calls
the same portion routine. At parameter λ on a grade-l curve, conformality
bounds each region's number of connected boundary pieces by a constant;
each piece has grade μ≤⌈βλ⌉ (tex 367). Lemma 2.3 bounds the region count
by O(2^(l−⌈βλ⌉)+1), so these resets cost

`O(2^l·λ·2^(−β²λ)) = O(2^l·λ·2^(−λ/25))`

per refinement round. The parameters strictly decrease (tex 379 and
`compute_visibility_map`), so summing over all rounds is bounded by
`O(2^l·sum_{λ>=1} λ·2^(−λ/25)) = O(2^l)`, hence O(n) for the whole curve.
Small grades and rounding contribute only constant factors.

Both series converge: for 0<r<1, `sum_{λ>=1} λ r^λ = r/(1−r)^2`.
The modified per-portion reset cost also remains below Lemma 4.1's
`O(λ² log λ·2^(74λ/75))` bound, since 1−β=4/5<74/75. Thus the reset
departure preserves **overall O(n)** in both phases, in the paper's
arithmetic-operation model. It still fails the local O(M) requirement.

## 4. LT79 Step 9: a parent-only ascent can miss the cycle

**Paper:** tex 626–628 directs an ascent from the triangle's third vertex y
through parent pointers until reaching the current cycle. The original
publication gives this instruction on p. 187 (PDF page 11). Step 8 permits
any non-tree edge and declares its heavier side inside (tex 603–614).

**Established contradiction:** the selected inside can contain the BFS root.
The root-to-y path then need not meet the fundamental cycle.
`tests/lt79-tests.cpp::test_root_inside_cycle` supplies a 13-vertex, 13-edge,
connected planar graph with explicit rotations. It has two faces of lengths
23 and 3, satisfying Euler's formula. The literal parent-only variant reaches
the root without encountering a cycle vertex and aborts. The current routine
passes both external separator and recursive-decomposition checks. This
fixture was reduced from a production §3.4 dual, where the same failure occurs.
The unrestricted Step 8 and parent-only Step 9 instructions are incompatible
for this valid input.

**Implemented correction:** one tree-path routine alternates ascents from y
and the cycle's minimum-depth vertex L. It either reaches the cycle directly,
or constructs y→lca(y,L)→L, including the required downward segment. There is
no failed-query recovery or substitute separator algorithm.

**Time:** if y's ascent hits the cycle after a steps, the unused ascent has
at most a+1 steps. Otherwise both searched legs belong to the constructed
path, with any overshoot bounded by the other leg's length. Search and cleanup
are O(path length+1). LT79 tex 661–663 charges each added path vertex to its
permanent removal from later interiors; tex 663–667 charges alternating edge
scans similarly. Thus Step 9 and the separator remain O(n). Assertions verify
that the constructed path consists of distinct off-cycle vertices connected
by tree edges and terminates on the cycle.

## 5. FM84 Algorithm 3: backtracking must skip the starting extremum

**Paper:** FM84 tex 432–435 saves `prev(current)`, removes `current`, then
chooses `next(first)` if the removed vertex equals `first`, or the saved
predecessor otherwise. The same instruction appears in the original publication
on p. 160 (PDF page 9), verified visually alongside the analysis on p. 161.
Tex 454–459 requires that neither height extremum be removed and
uses that invariant to prove correctness and linear time.

**Established contradiction:** take this clockwise simple polygon, with
`first = start = 0` and `last = 4`:

```
0: (0,8), 1: (6,6), 2: (4,4), 3: (1,2), 4: (0,0)
```

All heights are distinct, no adjacent edges are collinear, and the extrema
0 and 4 share the closing edge. The other chain strictly decreases in y,
so this meets the paper's unimonotone input definition. The initial current
vertex 1 is convex and emits triangle (0,1,2). Because 1 differs from `first`,
the printed backtrack selects its saved predecessor 0. With four vertices
remaining, vertex 0 is also convex, but its proposed triangle (4,0,2) strictly
contains vertex 3. Its three clockwise side determinants against vertex 3
are −8, −20, and −4. The two proposed triangles already have total doubled
area 48, exceeding the input's 44. The remaining three-vertex cycle has the
opposite winding. This is a failure on a valid input, not a transcription
error or an equal-height ambiguity.

`tests/triangulation-tests.cpp::check_endpoint_backtracking` checks the
contained vertex and requires the corrected triangle sequence
(0,1,2), (0,2,3), (0,3,4). The independent geometry checks verify exact area,
triangle interiors, and edge incidence.

**Implemented correction:** after deletion, select the saved predecessor
unless it is `start`; in that case select the new `next(start)` instead.
For clockwise input, `start` is the top extremum when the chain is on the
right and the bottom extremum when it is on the left, as tex 419–421 requires.
The current vertex stays strictly inside that chain. Assertions preserve
both extrema until the last triangle leaves a two-vertex list. There is no
alternative triangulation algorithm or failed-query recovery.

**Time:** for a piece of m vertices, there are m−2 deletions. Each deletion
causes at most one backward step. Counting current's position in the remaining
chain shows that every forward step is paired with a backward step, so there
are at most 2(m−2) convexity tests. Initializing and validating the list costs
O(m). Every triangle adds three adjacency entries, giving O(m) time and space.
Algorithm 2's total piece boundary length is at most 3n−6, so the combined
conversion and the complete algorithm retain overall O(n).
