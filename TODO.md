# TODO — Known Deferrals and Documented Deviations

Two kinds of entries, both cited against
`papers/chazelle1991-transcribed.tex`. **Deferrals**: work that depends
on machinery from later sections of the paper, each owned by a future
section. **Deviations**: places where the implementation deliberately
departs from the paper's literal text because the text is demonstrably
wrong; each carries the proof and the test that pins it. None is a
silent gap or an unexamined deviation.

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

## Blocked on future sections

1. **Merge wiring** (`merge.cpp::fuse`) — owned by **§3.4**.
   [C91 §3 tex 164]: "The merge proceeds in three stages" — stage 1 =
   fusion (§3.1, `fuse`), stage 2 = restoring conformality (§3.2,
   `restore_conformality`), stage 3 = maintaining granularity (§3.3,
   `enforce_granularity` + `Submap::normalize`); `merge()` runs exactly
   this pipeline. Stage 1's three-call wiring (pass 1, pass 2 with
   `junction_at_end = false`, `rebuild_submap`) is documented at the
   TODO in `fuse` but needs §3.4's real ray-shooting oracle before it
   can execute: `local_shoot` trips [C91 §3.1 tex 181]'s "local shoot
   must hit" assertion under stub oracles. This is now the ONLY
   blocker: stage 2 was implemented and wired 2026-07-03 and stage 3
   on 2026-07-04, so `merge()`'s contract ([C91 §3 tex 160]: "to merge
   S₁ and S₂ means to compute a normal-form γ-granular conformal
   submap of V(C)"; Lemma 3.5 tex 279) is met end-to-end once `fuse`
   produces a real fusion. Stages 2 and 3 gate themselves off while
   `fuse` still produces an empty submap.

2. **Wrapped (through-infinity) hits in the §3.1 fusion walk** — owned
   by **§3.4**. `RayHit.wrapped` ([C91 §2.1 tex 70]: a ray that misses
   everything "wraps around in the spherical plane until it hits C
   again") is fully supported by §3.2's `local_shoot_fused`
   (lexicographic (wrapped, distance) nearest order; d ≤ 0 behind the
   source), but §3.1's `local_shoot` and the case (i)/(ii) distance
   comparisons in `fuse_submaps` do not model the wrap metric yet —
   `local_shoot` asserts `!hit.wrapped` with a §3.4 TODO. A real fusion
   whose junction is a global y-extremum (a₀'s startup shot wraps —
   e.g. the apex outside-pair chord) needs this. Extend alongside
   §3.4's real oracle.

3. **End-to-end geometric validation of fusion** — owned by **§3.4**.
   The §3.1 fusion loop (both tour orientations) is validated
   line-by-line against the paper plus structural tests (sequence shape,
   startup region election, main-loop termination, dedup, full §2.4
   invariants on rebuild output), but has never run against a TRUE
   geometric visibility oracle. That residual risk is irreducible until
   §3.4 supplies one. (§3.2 now HAS run end-to-end against true
   geometric oracles — `tests/s3.2-tests.cpp`'s comb fixture — which is
   what exposed and fixed the inverted `shooting_direction` convention
   on 2026-07-03; see `src/merge/fusion.cpp`.)

4. **e2e `double_identify` completeness** — owned by **§4** (up-phase).
   The randomized test verifies soundness (every returned arc contains
   the query edge, counts ≤ 6 per [C91 §2.4 tex 144], non-empty at
   vertices) but not exhaustiveness against a symbolic reference; the
   conservative brute-force set is a superset by construction.
   Strengthening this requires a full symbolic V(C) reference map,
   natural once §4's up-phase provides one.
