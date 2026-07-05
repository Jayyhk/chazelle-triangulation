# TODO — Known Deferrals and Documented Deviations

Two kinds of entries, both cited against
`papers/chazelle1991-transcribed.tex`. **Deferrals**: work that depends
on machinery from later sections of the paper, each owned by a future
section. **Deviations**: places where the implementation deliberately
departs from the paper's literal text — either because the text is
demonstrably wrong, or because an equivalent representation better fits
the data structures; each carries the reasoning (and, where it fixes a
bug, the test that pins it). None is a silent gap or an unexamined
deviation.

## Documented paper deviations (not deferrals)

- **[C91 §3.3 tex 276]'s once-only claim** — the paper says "chords
  need be processed only once since the removals cannot make any chord
  removable if it was not already so before." That claim is provably
  half-wrong. Removability ([C91 §2.3 tex 121] criterion (ii),
  negated) is a conjunction: *incident upon a node of degree < 3* AND
  *contraction weight ≤ γ*.

  - **Weight half — the paper is right.** Contraction only glues and
    grows arcs, so a contraction weight that exceeds γ exceeds it
    forever; that half of removability is monotone and never needs
    re-checking.
  - **Degree half — the paper is wrong.** Contracting a leaf edge — a
    degree-1 node, which criterion (ii) explicitly covers ("degree
    less than 3") — leaves a merged node of degree d_u + 1 − 2 =
    d_u − 1. Degrees can DROP, so a chord examined while both
    endpoints had degree ≥ 3 can become removable later. No static
    processing order rescues the claim: parent-side degree drops can
    occur after any fixed position of a hub chord in the order.
  - **Why literal compliance breaks the paper's own theorem.** For
    γ ≥ total weight, γ-granularity ⟺ no chords at all (any surviving
    edge has a leaf endpoint: a degree-1 incidence with contraction
    weight ≤ γ, violating criterion (ii)). A single pass that visits a
    hub–hub chord before its leaf pockets leaves it behind, so the
    output is not γ-granular — contradicting Lemma 3.5 (tex 279),
    which everything downstream relies on (Lemma 2.3's O(n/γ + 1)
    region bound for merge outputs, the §4 up-phase budget). When the
    aside contradicts the theorem, the theorem wins.
  - **The fix (implemented in `enforce_granularity`).** Keep "process
    each exit chord in turn" verbatim; after each removal, re-examine
    only the ≤ 4 chords incident upon the merged node (only they can
    change status, and only via the degree drop). O(1) per removal, at
    most one removal per chord — tex 276's "time linear in the size of
    the submap tree" bound is preserved exactly, and the paper's
    "nondeterministic fashion" (any order, any valid removal set) is
    respected.
  - **Verification.** Mathematically as above, plus empirically
    (2026-07-04): with the re-check disabled, the criterion-(ii)
    postcondition assert inside `enforce_granularity` fires on the
    post-§3.2 six-peak comb with reversed chord order and γ = 100 — a
    removable hub chord survives the single pass; with the re-check,
    the fixed point is reached and all suites pass.
    `tests/s3.3-tests.cpp::test_degree_drop_recheck` pins the
    counterexample shape permanently.

- **[C91 §2.4 tex 142]'s double-backing arc-structure — we split
  instead.** An *equivalence*, not a correction: the paper's
  representation is correct; we choose a different-but-equivalent one so
  that the arc-sequence table is uniformly single-side.

  - **What the paper says.** A single arc that wraps around a C-endpoint
    ("double-backing") is stored as ONE arc-structure with two
    input-table pointers e₁, eₜ ([C91 §2.4 tex 142]: "an arc might wrap
    around both sides of C ... detected... as soon as we reach an edge
    of P incident upon an endpoint of C").
  - **What we do.** `rebuild_submap` (the §3.1 fusion-output builder)
    and `Submap::normalize` emit only single-side arc-structures
    (`first_side == last_side`): each wrap-spanning arc is split at the
    C-endpoint into a LEFT leg and a RIGHT leg — two separate table
    entries. (The `Arc` type can still express double-backing, and
    `double_identify` retains seam handling for it, but the merge/normal
    -form path never produces one.)
  - **Why.** The arc-sequence table is a total order — LEFT arcs
    ascending by `first_edge`, then RIGHT descending ([C91 §2.4(iii)
    tex 138]; enforced by `add_arc`). `double_identify` ([C91 §2.4
    tex 144]) breaks the circular sequence into those two linear halves
    and binary-searches each by (edge, y). A double-backing arc
    straddles both halves — it sorts into one run but is invisible to
    the other half's search — so the paper has to special-case the ≤ 2
    endpoint arcs at the seam. Splitting keeps every entry in exactly
    one half, so the O(log m) search — and every arc traversal, weight
    sum, and region-cycle walk — needs no per-arc "is this wrapped?"
    branch. (No recorded rationale from the original §2 author; this is
    the reconstructed engineering benefit, and it is the invariant the
    rest of the codebase is actually built on.)
  - **Cost / consequence.** ≤ 2 extra arc-structures per submap (one per
    C-endpoint wrap). Conformality still bounds a region to ≤ 4 PAPER
    arcs (degree ≤ 4, [C91 §2.3 tex 114]), but a region straddling both
    wraps then holds up to 4 + 2 = 6 arc-STRUCTURES.
    `fused_region_cycle` re-glues the legs into `LogicalArc`s wherever
    the paper's arc COUNT is meant (§3.2's "≤ 4 arcs" target);
    `is_conformal` (degree) is the paper invariant. Relied on by §2.4
    `double_identify`, §3.1 `collect_region_arcs` (the four C-endpoint
    arcs), §3.2 `LogicalArc`, and §3.4
    `RayShootingStructure::build_faces` (whose per-region assert is
    therefore ≤ 6, not ≤ 4 — the ≤ 4 form silently crashed on the first
    conformal fixture with a wrap-straddling region, 2026-07-04).
  - **Equivalence, not correction.** The paper itself conceptually
    splits the circular arc sequence into two linear sequences at the
    same two C-endpoints ([C91 §2.4 tex 144]); we make that split
    physical. Both representations encode the same V(C), and Lemma 2.3's
    bounds (O(n/γ) regions, O(γ) edges per region) are untouched by the
    O(1) extra structures.

## Blocked on future sections

1. **e2e `double_identify` completeness** — owned by **§4** (up-phase).
   The randomized test verifies soundness (every returned arc contains
   the query edge, counts ≤ 6 per [C91 §2.4 tex 144], non-empty at
   vertices) but not exhaustiveness against a symbolic reference; the
   conservative brute-force set is a superset by construction.
   Strengthening this requires a full symbolic V(C) reference map,
   natural once §4's up-phase provides one.

2. **Arc-cutting oracle** — owned by **§4** (up-phase).
   [C91 §3.4 tex 284]: "The arc-cutter is implemented by using the
   divide-and-conquer structure of the up-phase of the visibility
   algorithm ... we postpone the discussion of its implementation."
   Tests supply [C91 §3.0(ii) tex 170]-compliant cutters (validated by
   `assert_cut_postconditions`); the production cutter arrives with
   §4.1's chain decomposition ([C91 §4.1 tex 341–346]).

3. **[C91 §3.0(i) tex 169] oracle over piece structures** — owned by
   **§4** (up-phase). tex 169 specifies the report "in the absence of
   any obstacle except α'"; [C91 §4.1 tex 341] realizes that contract
   by decomposing α' into chains and shooting in each chain's own
   Lemma 3.6 structure (whose curve IS the piece). At the §3.4 stage
   the available structure is the one over Cᵢ itself
   (`SubmapRayShooter`), whose obstacle is Cᵢ ⊇ ᾱ'; the two reports
   agree at every §3.1/§3.2 call site (see the semantics note in
   `src/merge/ray_shooting.h`). When §4 lands, the merge oracles
   should be assembled per tex 341 and this note re-examined.
