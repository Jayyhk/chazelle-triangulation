# Documented Deviations from the Paper

Places where the implementation deliberately departs from the literal
text of `papers/chazelle1991-transcribed.tex` — either because the
text is demonstrably wrong (it contradicts the paper's own lemmas), or
because the paper's own conventions (SoS, the nonnull weight metric),
together with the implementation's null-length realization of
[C91 §4 tex 316]'s padding (a freedom the paper leaves open — tex 316
prescribes only "padding the curve with additional vertices, if
necessary", fixing neither placement nor edge lengths), force a
non-literal reading. Each entry carries
the tex-line citations, the reasoning, and — where the deviation fixes
a bug — the test that pins it. Nothing here is a silent gap or an
unexamined departure; anything NOT listed here matches the paper
exactly.

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
    endpoints had degree ≥ 3 can become removable later — and
    parent-side degree drops can occur after any fixed position of a
    hub chord in the order.  (A clairvoyant per-input order —
    processing a valid removal sequence first, after which criterion
    (ii) holds and every remaining chord is non-removable at its turn
    — would make one pass suffice; but computing such a sequence IS
    the fixpoint computation, and tex 276's "nondeterministic
    fashion" asserts that ANY order works, so a single bad order
    refutes the claim.)
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
  - **Counterexample (pinned).** On the post-§3.2 six-peak comb with
    reversed chord order and γ = 100, a single pass leaves a
    removable hub chord behind (the criterion-(ii) postcondition
    inside `enforce_granularity` fails without the re-check).
    `tests/s3.3-tests.cpp::test_degree_drop_recheck` pins the shape.

- **[C91 §3.2 tex 246]'s "which we ignore" for boundary subarcs** —
  the paper says of the ≤ 2 single-edge subarcs attached to A₁'s
  endpoints: "(which we ignore)".  Read as "never search them", that
  aside is provably wrong: when A₁ spans two partially-covered edges,
  the arc-cutter's decomposition consists ONLY of boundary pieces, and
  Lemma 3.3's guaranteed vertex ([C91 §3.2 tex 256–264]) lives there —
  skipping them would strand the tex-264 chord-addition loop on a
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

- **[C91 §4.1 tex 341]'s piece sizes under padding** — tex 341 says
  each vertex-to-vertex piece of the arc-cutter's decomposition has
  "at most 2^{⌈βλ⌉} edges", and keeps every cover chain "in
  grades at most ⌈βλ⌉" (same line).  Read against the
  implementation's null-length realization of [C91 §4 tex 316]'s
  padding (chain.cpp pad_curve appends coincident vertices as an
  all-null suffix; tex 316 itself prescribes neither placement nor
  edge lengths) and [C91 §2.2 tex 106] ("weights count only edges of
  nonzero length"), the edge counts are NONNULL counts: a piece with
  ≤ 2^{⌈βλ⌉} nonnull edges may span arbitrarily many null padding
  edges, so its dyadic cover admits ALL-NULL chains of any processed
  grade < λ.  Under that realization this is forced, not merely
  convenient: the premise itself ("any arc α ... consists of at most
  γ edges", tex 341) is DERIVED from γ-granularity, whose weight
  metric counts nonnull edges only — with null padding edges present
  the literal all-edges reading of tex 341 is unsatisfiable — and
  covering a null run with grade ≤ ⌈βλ⌉ chains would need up to
  2^{λ−⌈βλ⌉} of them, destroying tex 341's g(γ) = O(λ).  The h-bound
  survives by the paper's own default clause ([C91 §2.3 tex 123]:
  "if (i) holds but
  the submap has no exit chord, it is still said to be γ-granular"):
  an all-null chain's canonical submap is chordless with total
  weight 0, hence γ-granular for EVERY γ, and its single-leaf tree
  decomposition keeps it searchable by Lemma 3.2's descent.  Every
  chain containing a nonnull edge stays in grade ≤ ⌈βλ⌉ exactly as
  tex 341 states.  Implemented in visibility/up_phase.cpp
  (piece_cover).

- **All-null dyadic parts attach by extension, not merge** — Lemma
  4.1 (tex 336, proof 339–341) partitions a portion v_a..v_b into ≤ 2λ chains
  and "merges them two-by-two".  A suffix part consisting ONLY of
  null padding edges ([C91 §4 tex 316]) contributes no visibility
  information usable by a canonical submap: the full maps DIFFER
  (the tail vertices carry their own null-cluster chords in V), but
  a submap may omit any subset of V's chords, and a canonical
  submap of V(D ∪ tail) is obtained from one of V(D) directly — the
  tail's chords are omitted, every kept chord remains a visible
  pair (the tail occupies a single point, so only sightlines
  through the old boundary point are affected — exactly the
  relabeled incidences), weights are unchanged (null edges,
  [C91 §2.2 tex 106], with arc.h's zero-coverage counting at the
  old corner), and conformality/granularity carry over (asserted).
  Merging such a part through §3.1's fusion instead would put an
  all-null curve on one side of the junction, a maximally
  degenerate configuration (every raw distance 0) the walk never
  needs elsewhere.  The implementation therefore attaches all-null
  parts by EXTENSION (visibility/up_phase.cpp,
  extend_with_null_tail): move the submap's end vertex across the
  null tail, relabel the bounded set of chords/arcs anchored at the
  old boundary when it becomes a y-extremum of the extended curve,
  then refresh edge counts and the tree decomposition.  Cost:
  O(M log M) per extension, M = the submap size (the chord/arc
  relabel and count refresh are O(M); rebuilding the tree
  decomposition is O(M log M), [C91 §2.3 tex 116]) — a λ factor over
  Lemma 4.1's granularity-reset step O(2^{λ(1−β)}) (tex 339) and
  still dominated by the lemma's own merge bound (tex 348–350), so
  the lemma's time bound is unaffected.  No merge ever sees an all-null input
  or a junction at a padding vertex (all-null parts form a suffix
  of the dyadic part sequence, so merge junctions sit at boundaries
  between real-containing parts, i.e. at real vertices).
- **[C91 §3.1 tex 206]'s case (i) action vs invariant (B)** —
  tex 206's case (i) action says "set p = a_j, let the current region
  still be R".  Taken literally, that action can break the paper's
  own loop invariant: when the case (i) hit s lands EXACTLY on an
  exit-chord endpoint of S₂ at the sightline's level, the new p's
  sight rides that chord, and invariant (B) (tex 195) prescribes the
  tie rule "if p lies on a chord between two regions of S₂, the
  current region should be the one that we enter as we locally leave
  p clockwise around ∂C₁" — which may be the OTHER side of the chord.
  Keeping R unchanged then strands the walk on the wrong side, and
  the Lemma 3.1 trichotomy fails at a later stop (its correctness
  proof's "the point a' ... cannot be equal to p" presumes the
  election), leaving stale wrap chords un-revalidated (tex 224).
  The implementation applies the invariant-B election in the case (i)
  action exactly as the startup does (tex 191; fusion.cpp, the
  "invariant-B tie" block); the endpoint match is POSITIONAL
  (crossing vertex + side, not exact edge label — [C91 §3.0(i)
  tex 169]: a point at a shared crossing vertex is contained in both
  incident edges, so the recorded label may be the neighbor edge's).  When the summary action contradicts the
  loop invariant, the invariant wins.  Pinned by
  tests/s4.1-tests.cpp::test_deep_grades_invariant_b (a 248-vertex
  curve whose [208,224] merge hits the configuration).
