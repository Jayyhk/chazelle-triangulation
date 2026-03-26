import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
import os
from matplotlib.path import Path


def intersect(p1, p2, p3, p4):
    """Returns True if line segment p1p2 and p3p4 intersect (excluding endpoints)."""

    def ccw(A, B, C):
        return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])

    # Standard segment intersection check
    # Check if segments intersect at all
    if ccw(p1, p3, p4) != ccw(p2, p3, p4) and ccw(p1, p2, p3) != ccw(p1, p2, p4):
        # To exclude endpoints, we check if any endpoint lies on the other segment
        # However, for triangulation, we specifically care about INTERIOR intersections.
        return True
    return False


def edge_x_at_y(verts, edge_idx, y):
    """Return x coordinate where polygon edge (edge_idx -> edge_idx+1) crosses y."""
    n = len(verts)
    v1 = verts[edge_idx]
    v2 = verts[(edge_idx + 1) % n]
    y1, y2 = v1["y"], v2["y"]
    if abs(y2 - y1) < 1e-12:
        return (v1["x"] + v2["x"]) / 2
    t = (y - y1) / (y2 - y1)
    return v1["x"] + t * (v2["x"] - v1["x"])


def is_edge_illegal(p1, p2, poly_pts, edge_indices, all_verts):
    # 1. Check midpoint
    mid = [(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2]
    if not Path(poly_pts).contains_point(mid):
        return True

    # 2. Check for intersections with ANY boundary edge
    # Boundary edges are (i, i+1)
    n = len(poly_pts)
    for i in range(n):
        v1_idx = i
        v2_idx = (i + 1) % n

        # Skip if the boundary edge shares a vertex with the test edge
        if v1_idx in edge_indices or v2_idx in edge_indices:
            continue

        b1 = poly_pts[v1_idx]
        b2 = poly_pts[v2_idx]

        if intersect(p1, p2, b1, b2):
            return True

    return False


def _plot_common(ax, verts, poly_pts):
    """Plot polygon boundary, vertices, and indices."""
    boundary = patches.Polygon(
        poly_pts, closed=True, fill=False, edgecolor="black", linewidth=2.5, zorder=30
    )
    ax.add_patch(boundary)
    xs = [v["x"] for v in verts]
    ys = [v["y"] for v in verts]
    ax.scatter(xs, ys, color="black", s=20, zorder=35)
    for i, v in enumerate(verts):
        ax.annotate(
            str(v.get("i", i)),
            (v["x"], v["y"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
            fontweight="bold",
            color="blue",
            zorder=40,
            bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.6, ec="none"),
        )
    ax.set_aspect("equal")
    plt.grid(True, linestyle="--", alpha=0.3)


def visualize_triangulation(data, output_path):
    """Produce triangulation comparison image (reference + Chazelle)."""
    verts = data["vertices"]
    ref_tris = data.get("reference_triangulation", [])
    chaz_tris = data.get("chazelle_triangulation", [])
    name = data.get("name", "unknown")
    poly_pts = [[v["x"], v["y"]] for v in verts]

    fig, ax = plt.subplots(figsize=(14, 14))

    # Plot Reference Triangles (Very Light Green Fill, No Edges)
    for t in ref_tris:
        p1 = [verts[t[0]]["x"], verts[t[0]]["y"]]
        p2 = [verts[t[1]]["x"], verts[t[1]]["y"]]
        p3 = [verts[t[2]]["x"], verts[t[2]]["y"]]
        poly = patches.Polygon(
            [p1, p2, p3],
            closed=True,
            facecolor="tab:green",
            alpha=0.1,
            fill=True,
            edgecolor="none",
            zorder=5,
        )
        ax.add_patch(poly)

    # Plot Chazelle Triangles with Robust Edge Classification
    seen_edges = set()
    for t in chaz_tris:
        edges = [
            tuple(sorted((t[0], t[1]))),
            tuple(sorted((t[1], t[2]))),
            tuple(sorted((t[2], t[0]))),
        ]
        for e in edges:
            if e in seen_edges:
                continue
            seen_edges.add(e)
            p1 = [verts[e[0]]["x"], verts[e[0]]["y"]]
            p2 = [verts[e[1]]["x"], verts[e[1]]["y"]]
            is_boundary = any(
                tuple(sorted((i, (i + 1) % len(verts)))) == e for i in range(len(verts))
            )
            if is_boundary:
                continue
            if is_edge_illegal(p1, p2, poly_pts, e, verts):
                color, lw, alpha, z = "tab:red", 2.0, 1.0, 12
            else:
                color, lw, alpha, z = "#d1d5db", 1.0, 0.6, 12
            ax.plot(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                color=color,
                linewidth=lw,
                alpha=alpha,
                zorder=z,
            )

    _plot_common(ax, verts, poly_pts)
    ax.set_title(
        f"Triangulation: {name}\nGrey = Chazelle Edges, Red = Illegal", fontsize=14
    )
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path}")


def visualize_trapezoids(data, output_path):
    """Produce trapezoidalization comparison image."""
    verts = data["vertices"]
    chaz_traps = data.get(
        "chazelle_trapezoidalation", data.get("trapezoidalization", [])
    )
    ref_traps = data.get("reference_trapezoidalation", [])
    name = data.get("name", "unknown")
    poly_pts = [[v["x"], v["y"]] for v in verts]

    fig, ax = plt.subplots(figsize=(14, 14))

    def traps_equal(t1, t2):
        return (
            t1["top_vertex"] == t2["top_vertex"]
            and t1["bottom_vertex"] == t2["bottom_vertex"]
            and t1["left_edge"] == t2["left_edge"]
            and t1["right_edge"] == t2["right_edge"]
        )

    def draw_traps(traps, ref_traps, base_z=50):
        for t in traps:
            # Check if this Chazelle trapezoid exists in the reference set
            is_correct = any(traps_equal(t, rt) for rt in ref_traps)

            color = "#00FF00" if is_correct else "tab:red"
            alpha = 1.0 if is_correct else 0.8
            linestyle = "-"

            top_v = t["top_vertex"]
            bot_v = t["bottom_vertex"]
            le = t["left_edge"]
            re = t["right_edge"]
            y_top = verts[top_v]["y"]
            y_bot = verts[bot_v]["y"]
            x_lt = edge_x_at_y(verts, le, y_top)
            x_rt = edge_x_at_y(verts, re, y_top)
            x_lb = edge_x_at_y(verts, le, y_bot)
            x_rb = edge_x_at_y(verts, re, y_bot)

            # Top/Bottom horizontal lines
            ax.plot(
                [x_lt, x_rt],
                [y_top, y_top],
                color=color,
                linestyle=linestyle,
                alpha=alpha,
                zorder=base_z,
            )
            ax.plot(
                [x_lb, x_rb],
                [y_bot, y_bot],
                color=color,
                linestyle=linestyle,
                alpha=alpha,
                zorder=base_z,
            )

            # Left/Right walls
            ax.plot(
                [x_lt, x_lb],
                [y_top, y_bot],
                color=color,
                linestyle=linestyle,
                alpha=alpha,
                zorder=base_z,
            )
            ax.plot(
                [x_rt, x_rb],
                [y_top, y_bot],
                color=color,
                linestyle=linestyle,
                alpha=alpha,
                zorder=base_z,
            )

    # Draw Chazelle traps only, color-coded by correctness
    draw_traps(chaz_traps, ref_traps, base_z=50)

    _plot_common(ax, verts, poly_pts)
    ax.set_title(
        f"Trapezoidalization: {name}\nBright Green = Correct, Red = Incorrect",
        fontsize=14,
    )
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path}")


def _sanitize_name(name):
    """Make name safe for use in filenames."""
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in str(name))


def visualize(
    json_path, output_dir=None, use_name_in_output=False, trapezoid_only=False
):
    """Load JSON and produce triangulation and/or trapezoid comparison images."""
    if not os.path.exists(json_path):
        print(f"Error: {json_path} not found.")
        return

    with open(json_path, "r") as f:
        data = json.load(f)

    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(json_path))
    if output_dir == "":
        output_dir = os.path.dirname(os.path.abspath(__file__))

    if use_name_in_output:
        safe_name = _sanitize_name(data.get("name", "unknown"))
        tri_path = os.path.join(output_dir, f"triangulation_comparison_{safe_name}.png")
        trap_path = os.path.join(output_dir, f"trapezoid_comparison_{safe_name}.png")
    else:
        tri_path = os.path.join(output_dir, "triangulation_comparison.png")
        trap_path = os.path.join(output_dir, "trapezoid_comparison.png")

    if not trapezoid_only:
        visualize_triangulation(data, tri_path)
    visualize_trapezoids(data, trap_path)


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, "failed_polygon.json")
    use_name = "--named" in sys.argv
    trapezoid_only = "--trapezoid-only" in sys.argv
    for arg in sys.argv[1:]:
        if not arg.startswith("-"):
            path = arg
            break
    visualize(path, use_name_in_output=use_name, trapezoid_only=trapezoid_only)
