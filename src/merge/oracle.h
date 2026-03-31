#pragma once

/// [C91 §3]: Oracle primitives for submap merging.
///
/// "To facilitate the exposition we assume that we have at our
/// disposal two primitives."
///
/// These are abstract interfaces — implementations are provided
/// by §3.4 and the up-phase (§4).

#include "../polygon/point.h"
#include "../polygon/polygon.h"
#include "../submap/submap.h"

#include <cstddef>
#include <vector>

namespace chazelle {

/// [C91 §3]: A subarc specification — a contiguous piece of an arc
/// on ∂C, identified by its two endpoint edges/sides in the input
/// table.  Used as input to both oracles.
struct Subarc {
    std::size_t first_edge;
    Side first_side;
    std::size_t last_edge;
    Side last_side;
};

/// [C91 §3]: Result of a ray-shooting query.
struct RayHit {
    bool hit = false;          ///< True if the ray hit the subarc.
    double x = 0.0;            ///< x-coordinate of the hit point.
    double y = 0.0;            ///< y-coordinate of the hit point.
    std::size_t edge = 0;      ///< Edge of P containing the hit point.
    Side side = LEFT;          ///< Which side of ∂C.
};

/// [C91 §3]: Result of an arc-cutting query — one piece.
struct ArcPiece {
    Subarc subarc;             ///< The piece's location on ∂C.
    const Submap* submap;      ///< Conformal submap of V(piece), or
                               ///< nullptr for single-edge pieces.
    const Polygon* curve;      ///< The piece's underlying curve.
};

/// [C91 §3]: "a ray-shooting oracle, which allows us to shoot a
/// horizontal ray toward any subarc of S₁ or S₂ and discover which
/// point, if any, is first hit by the ray; this gives us a way to
/// discover new chords."
///
/// Detailed specification: §3.0(i).
/// Implementation: §3.4.
struct RayShootingOracle {
    virtual ~RayShootingOracle() = default;

    /// Shoot a horizontal ray from point p in the given direction
    /// toward subarc α'.  Returns the first point hit (if any).
    virtual RayHit shoot(Point p, Side direction,
                         const Subarc& target) const = 0;
};

/// [C91 §3]: "an arc-cutting oracle, which enables us to cut up any
/// subarc into a small number of pieces for which conformal submaps
/// of the appropriate granularity have already been computed.  This
/// is to be used for restoring the conformality of merged submaps."
///
/// Detailed specification: §3.0(ii).
/// Implementation: §4 (up-phase).
struct ArcCuttingOracle {
    virtual ~ArcCuttingOracle() = default;

    /// Subdivide subarc α' into pieces, each with a precomputed
    /// conformal submap of the appropriate granularity.
    virtual std::vector<ArcPiece> cut(const Subarc& target) const = 0;
};

} // namespace chazelle
