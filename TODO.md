# TODO — Known Deferrals

Items examined during the 2026-07-03 strict audits (§2 and §3.0–§3.1) and
deliberately deferred or left as-is. Each is documented, cited against
`papers/chazelle1991-transcribed.tex`, and owned by a future section.
None is a silent gap or an unexamined deviation.

## Blocked on future sections

1. **Merge wiring** (`merge.cpp::fuse`) — owned by **§3.4**. The
   three-call wiring (pass 1, pass 2 with `junction_at_end = false`,
   `rebuild_submap`) is documented at the TODO but needs §3.4's real
   ray-shooting oracle before it can execute: `local_shoot` trips
   [C91 §3.1 tex 181]'s "local shoot must hit" assertion under stub
   oracles, and `merge()`'s contract (normal-form γ-granular conformal
   submap) is unmeetable while §3.2/§3.3 are stubs anyway. The former
   blockers (b) symmetric-pass orientation and (c) inventory dedup were
   implemented 2026-07-03.

2. **End-to-end geometric validation of fusion** — owned by **§3.4**.
   The §3.1 fusion loop (both tour orientations) is validated
   line-by-line against the paper plus structural tests (sequence shape,
   startup region election, main-loop termination, dedup, full §2.4
   invariants on rebuild output), but has never run against a TRUE
   geometric visibility oracle. That residual risk is irreducible until
   §3.4 supplies one.

3. **Orphaned arcs after cascaded removals** — owned by **§3.3's
   `normalize()`** (tex 276 *"We can now put S in normal form"*). Arcs
   whose bounding chords are all removed drop out of every chord's
   adjacency list, so `region_weight` can under-read until the submap is
   re-normalized. The s2 e2e suite documents each skipped check with a
   TODO. Not reachable in currently implemented paper flows (fusion
   output carries chords at C's endpoints; §3.3 is not yet on disk).

4. **e2e `double_identify` completeness** — owned by **§4** (up-phase).
   The randomized test verifies soundness (every returned arc contains
   the query edge, counts ≤ 6 per [C91 §2.4 tex 144], non-empty at
   vertices) but not exhaustiveness against a symbolic reference; the
   conservative brute-force set is a superset by construction.
   Strengthening this requires a full symbolic V(C) reference map,
   natural once §4's up-phase provides one.

## Examined and deliberately left as-is

5. **Same-position enumeration ties** (outside-pair duplicates at an
   extremum, [C91 §2.1 tex 72]): `build_fusion_sequence`'s comparator
   cannot order two stops at the identical ∂C point; their relative
   order does not affect any case test (A_j between them is zero-length)
   — revisit if §3.4/§4 machinery becomes sensitive to intra-point order.

6. **Boundary-edge count convention** at mid-edge splits (both arcs count
   the shared split edge, [C91 §2.4 tex 133]) — documented §2-audit
   convention, cache-consistent on both sides of every weight
   comparison; vertex-position endpoints no longer participate (fixed in
   the §3.0–3.1 audit's rebuild-sweep cursor handling).

7. **`s_dist == t_dist` tie in case (i)** ([C91 §3.1 tex 220]) — strict
   `<` is correct: equality would require the C₁- and C₂-hits to
   coincide, possible only at the junction vertex, where the old S₁
   chord remains valid in V(C) and carries over (tex 224's junction
   category). No action needed; recorded so the reasoning isn't redone.
