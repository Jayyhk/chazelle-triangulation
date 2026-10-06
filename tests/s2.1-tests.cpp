#include "algorithm/polygon/perturbation.h"
#include "algorithm/polygon/polygon.h"

#include <cassert>
#include <cstdio>

using namespace chazelle;

static void test_is_endpoint() {
    Polygon tri({{0.0, 0.0, 0}, {4.0, 1.0, 1}, {2.0, 3.0, 2}});
    assert(tri.is_endpoint(0));
    assert(!tri.is_endpoint(1));
    assert(tri.is_endpoint(2));

    Polygon curve({{0, 0, 0}, {1, 2, 1}, {3, 1, 2}, {5, 4, 3}, {6, 0, 4}});
    assert(curve.is_endpoint(0));
    assert(!curve.is_endpoint(1));
    assert(!curve.is_endpoint(2));
    assert(!curve.is_endpoint(3));
    assert(curve.is_endpoint(4));

    std::printf("  [PASS] is_endpoint\n");
}

static void test_local_extremum_fns() {
    Point a{0, 5, 0}, b{2, 1, 1}, c{4, 5, 2};
    assert(is_local_y_minimum(a, b, c));
    assert(!is_local_y_maximum(a, b, c));
    assert(is_local_y_extremum(a, b, c));

    Point d{0, 1, 0}, e{2, 5, 1}, f{4, 1, 2};
    assert(!is_local_y_minimum(d, e, f));
    assert(is_local_y_maximum(d, e, f));
    assert(is_local_y_extremum(d, e, f));

    Point g{0, 0, 0}, h{1, 1, 1}, i{2, 2, 2};
    assert(!is_local_y_minimum(g, h, i));
    assert(!is_local_y_maximum(g, h, i));
    assert(!is_local_y_extremum(g, h, i));

    Point j{0, 2, 0}, k{1, 1, 1}, l{2, 0, 2};
    assert(!is_local_y_minimum(j, k, l));
    assert(!is_local_y_maximum(j, k, l));
    assert(!is_local_y_extremum(j, k, l));

    std::printf("  [PASS] local_extremum_fns\n");
}

static void test_sos_extremum() {
    Point a{0, 0, 0}, b{1, 0, 1}, c{2, 0, 2};
    assert(!is_local_y_extremum(a, b, c));

    assert(is_local_y_maximum(c, a, b));

    assert(is_local_y_minimum(a, c, b));

    Point d{0, 5, 0}, e{1, 5, 1}, f{2, 0, 2};
    assert(!is_local_y_extremum(d, e, f));

    assert(is_local_y_maximum(e, d, f));

    std::printf("  [PASS] sos_extremum\n");
}

static void test_polygon_is_y_extremum() {
    Polygon v_shape({{0, 0, 0}, {2, 3, 1}, {4, 0, 2}});
    assert(!v_shape.is_y_extremum(0));
    assert(v_shape.is_y_extremum(1));
    assert(!v_shape.is_y_extremum(2));

    Polygon w_shape({{0, 4, 0}, {1, 0, 1}, {2, 3, 2}, {3, 0, 3}, {4, 4, 4}});
    assert(!w_shape.is_y_extremum(0));
    assert(w_shape.is_y_extremum(1));
    assert(w_shape.is_y_extremum(2));
    assert(w_shape.is_y_extremum(3));
    assert(!w_shape.is_y_extremum(4));

    Polygon mono({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}});
    assert(!mono.is_y_extremum(0));
    assert(!mono.is_y_extremum(1));
    assert(!mono.is_y_extremum(2));
    assert(!mono.is_y_extremum(3));

    std::printf("  [PASS] polygon_is_y_extremum\n");
}

static void test_fig22_cases() {
    Polygon curve({{0, 0, 0}, {1, 2, 1}, {2, 5, 2}, {3, 3, 3}, {4, 1, 4}});

    assert(curve.is_endpoint(0) && !curve.is_y_extremum(0));
    assert(curve.is_endpoint(4) && !curve.is_y_extremum(4));

    assert(!curve.is_endpoint(2) && curve.is_y_extremum(2));

    assert(!curve.is_endpoint(1) && !curve.is_y_extremum(1));
    assert(!curve.is_endpoint(3) && !curve.is_y_extremum(3));

    std::printf("  [PASS] fig22_cases\n");
}

static void test_side_enum() {
    assert(LEFT == 0);
    assert(RIGHT == 1);
    assert(LEFT != RIGHT);

    Side s = LEFT;
    assert(s == LEFT);
    s = RIGHT;
    assert(s == RIGHT);

    std::printf("  [PASS] side_enum\n");
}

static void test_polygon_construction() {
    Polygon p({{0, 0, 0}, {4, 0, 1}, {4, 4, 2}, {0, 4, 3}});
    assert(p.num_vertices() == 4);
    assert(p.num_edges() == 3);

    assert(p.edge(0).start_idx == 0 && p.edge(0).end_idx == 1);
    assert(p.edge(1).start_idx == 1 && p.edge(1).end_idx == 2);
    assert(p.edge(2).start_idx == 2 && p.edge(2).end_idx == 3);

    assert(p.vertex(0).x == 0.0 && p.vertex(0).y == 0.0);
    assert(p.vertex(2).x == 4.0 && p.vertex(2).y == 4.0);

    std::printf("  [PASS] polygon_construction\n");
}

static void test_endpoint_crossing_coordinates() {
    for (bool descending : {false, true}) {
        Polygon c({{32.0, descending ? 1.0 : 0.0, 0}, {1.1, descending ? 0.0 : 1.0, 1}});
        for (std::size_t v : {std::size_t{0}, std::size_t{1}}) {
            const SymbolicY sy = symbolic_y_of(c.vertex(v));
            Exact x = 0.0;
            assert(edge_x_at_y(c, 0, sy) == c.vertex(v).x);
            assert(edge_crossing_x(c, 0, sy, &x) && x == c.vertex(v).x);
        }
        const SymbolicY tied{c.vertex(1).y, descending ? std::size_t{0} : std::size_t{2}};
        Exact x = 0.0;
        assert(edge_x_at_y(c, 0, tied) == c.vertex(1).x);
        assert(edge_crossing_x(c, 0, tied, &x) && x == c.vertex(1).x);
    }
    std::printf("  [PASS] endpoint_crossing_coordinates\n");
}

int main() {
    std::printf("[C91 §2.1 tests]:\n");
    test_is_endpoint();
    test_local_extremum_fns();
    test_sos_extremum();
    test_polygon_is_y_extremum();
    test_fig22_cases();
    test_side_enum();
    test_polygon_construction();
    test_endpoint_crossing_coordinates();
    std::printf("All §2.1 tests passed.\n");
    return 0;
}
