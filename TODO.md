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

- **[C91 §3.2 tex 248]'s "which we ignore" for boundary subarcs** —
  the paper says of the ≤ 2 single-edge subarcs attached to A₁'s
  endpoints: "(which we ignore)".  Read as "never search them", that
  aside is provably wrong: when A₁ spans two partially-covered edges,
  the arc-cutter's decomposition consists ONLY of boundary pieces, and
  Lemma 3.3's guaranteed vertex ([C91 §3.2 tex 256–264]) lives there —
  skipping them would strand the tex-267 chord-addition loop on a
  region it can never fix, contradicting Lemma 3.4.  The
  implementation follows the reading under which the paper's own
  lemmas hold: the parenthetical scopes the sentence it sits in (the
  boundary pieces have no submaps, so they are outside the
  tree-descent machinery), while the next sentence's "we search each
  subarc in turn" includes them — their O(1) on-A₁ vertices are
  checked directly (conformality.cpp), costing O(1) oracle shots per
  piece, within Lemma 3.2's O(f(γ₂)g(γ₂)(h(γ₂) + log γ₂)) budget.
  When the aside contradicts the theorem, the theorem wins.  Pinned by
  tests/s3.2-tests.cpp test 6(a).

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

3. **`endpoint_shot`'s non-wrapping-chord assumption** — owned by
   **§4** (up-phase).  [C91 §3.2 tex 246] descends S_α's tree
   decomposition testing "a chord ab of S_α" with no restriction, and
   an h-granular conformal S_α CAN retain chords that run through
   infinity ([C91 §2.1 tex 70] — e.g. a global-extremum outside cap
   whose contraction would over-weight the merged arc survives
   criterion (ii)).  `endpoint_shot` (conformality.cpp) computes
   a′/b′ with plain signed-x distances and asserts a direct chord
   (`(dir == RIGHT) == (other.x > from.x)`, `!shot.hit.wrapped`).
   Unreachable today: every S_α on main comes from the §3.0(ii) test
   cutters, which supply chordless (or single-direct-chord) submaps.
   When §4.1's production arc-cutter lands, the a′/b′ arithmetic must
   move to the lexicographic (wrapped, distance) traversal metric used
   by local_shoot_fused/fusion, and the two asserts must admit
   wrapping ab (Lemma 2.4, tex 150–156, is purely topological and
   already covers it).

4. **[C91 §3.0(i) tex 169] oracle over piece structures** — owned by
   **§4** (up-phase). tex 169 specifies the report "in the absence of
   any obstacle except α'"; [C91 §4.1 tex 341] realizes that contract
   by decomposing α' into chains and shooting in each chain's own
   Lemma 3.6 structure (whose curve IS the piece). At the §3.4 stage
   the available structure is the one over Cᵢ itself
   (`SubmapRayShooter`), whose obstacle is Cᵢ ⊇ ᾱ'; the two reports
   agree at every §3.1/§3.2 call site (see the semantics note in
   `src/merge/ray_shooting.h`). When §4 lands, the merge oracles
   should be assembled per tex 341 and this note re-examined.
