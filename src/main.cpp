#include "geometry/point.h"
#include "geometry/polygon.h"
#include "chazelle/grade_storage.h"
#include "chazelle/up_phase.h"
#include "chazelle/down_phase.h"
#include "triangulate/convert.h"
#include "triangulate/trapezoid.h"
#include "triangulate/vertex_node.h"
#include "triangulate/triangulate.h"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <vector>

int main() {
    // ── Read polygon from stdin ─────────────────────────────────
    // Format: first line = vertex count n,
    //         then n lines of "x y" coordinates in boundary order.
    std::size_t n = 0;
    if (!(std::cin >> n) || n < 3) {
        std::cerr << "error: need at least 3 vertices\n";
        return EXIT_FAILURE;
    }

    std::vector<chazelle::Point> vertices(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (!(std::cin >> vertices[i].x >> vertices[i].y)) {
            std::cerr << "error: failed to read vertex " << i << '\n';
            return EXIT_FAILURE;
        }
    }

    // ── Build polygon (pads to 2^p + 1 vertices) ───────────────
    chazelle::Polygon polygon(std::move(vertices));

    std::cout << "polygon: " << polygon.original_size()
              << " vertices (padded to " << polygon.num_vertices()
              << "), " << polygon.num_edges()
              << " edges, " << polygon.num_grades()
              << " grades\n";

    auto t0 = std::chrono::steady_clock::now();

    // ── Up-phase ────────────────────────────────────────────────
    chazelle::GradeStorage storage;
    chazelle::up_phase(polygon, storage);

    auto t1 = std::chrono::steady_clock::now();
    std::cout << "up-phase:   "
              << std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()
              << " µs\n";

    // ── Down-phase ──────────────────────────────────────────────
    chazelle::Submap vp = chazelle::down_phase(polygon, storage);

    auto t2 = std::chrono::steady_clock::now();
    std::cout << "down-phase: "
              << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
              << " µs\n";

    std::cout << "V(P): " << vp.num_nodes() << " regions, "
              << vp.num_chords() << " chords, "
              << vp.num_arcs() << " arcs\n";

    // ── Convert V(P) → Fournier-Montuno format ─────────────────
    // Each region of V(P) becomes a trapezoid; each polygon vertex
    // becomes a node in the vertex-node linked list.
    auto fm = chazelle::convert_submap_to_fm(vp, polygon);

    // ── Fournier-Montuno triangulation ──────────────────────────
    auto triangles = chazelle::fm_triangulate(
        fm.nodes, fm.traps, fm.first, fm.last);

    auto t3 = std::chrono::steady_clock::now();
    std::cout << "triangulate: "
              << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count()
              << " µs\n";

    // ── Output ──────────────────────────────────────────────────
    std::cout << triangles.size() << " triangles\n";
    for (auto& [a, b, c] : triangles) {
        std::cout << a << ' ' << b << ' ' << c << '\n';
    }

    auto total = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t0).count();
    std::cout << "total: " << total << " µs\n";

    return EXIT_SUCCESS;
}
