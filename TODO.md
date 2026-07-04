# TODO — Known Deferrals

Work deferred during because it depends on machinery from later sections 
of the paper. Each is documented, cited against 
`papers/chazelle1991-transcribed.tex`, and owned by a future section. 
None is a silent gap or an unexamined deviation.

## Blocked on future sections

1. **Merge wiring** (`merge.cpp::fuse`) — owned by **§3.4**. The
   three-call wiring (pass 1, pass 2 with `junction_at_end = false`,
   `rebuild_submap`) is documented at the TODO but needs §3.4's real
   ray-shooting oracle before it can execute: `local_shoot` trips
   [C91 §3.1 tex 181]'s "local shoot must hit" assertion under stub
   oracles, and `merge()`'s contract (normal-form γ-granular conformal
   submap) is unmeetable while §3.3 is a stub anyway. The former
   blockers (b) symmetric-pass orientation and (c) inventory dedup were
   implemented 2026-07-03; the §3.2 stage (`restore_conformality`,
   `src/merge/conformality.{h,cpp}`) was implemented and wired
   2026-07-03 (the merge pipeline gates it off while `fuse` still
   produces an empty submap).

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

4. **Orphaned arcs after cascaded removals** — owned by **§3.3's
   `normalize()`** (tex 276 *"We can now put S in normal form"*). Arcs
   whose bounding chords are all removed drop out of every chord's
   adjacency list, so `region_weight` can under-read until the submap is
   re-normalized. The s2 e2e suite documents each skipped check with a
   TODO. Not reachable in currently implemented paper flows (fusion
   output carries chords at C's endpoints; §3.3 is not yet on disk).

5. **e2e `double_identify` completeness** — owned by **§4** (up-phase).
   The randomized test verifies soundness (every returned arc contains
   the query edge, counts ≤ 6 per [C91 §2.4 tex 144], non-empty at
   vertices) but not exhaustiveness against a symbolic reference; the
   conservative brute-force set is a superset by construction.
   Strengthening this requires a full symbolic V(C) reference map,
   natural once §4's up-phase provides one.

6. **Normal form after §3.2 insertions** — owned by **§3.3** (tex 276).
   `Submap::insert_chord` appends split-arc halves at the table's end,
   breaking canonical arc-sequence order ([C91 §2.4(iii) tex 138]);
   `compacted_` is cleared so `double_identify` fails fast. §3.3's
   normalize step re-sorts the table and rebuilds the tree
   decomposition, exactly as the paper re-establishes normal form after
   conformality restoration.
