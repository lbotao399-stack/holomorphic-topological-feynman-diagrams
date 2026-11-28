#!/usr/bin/env python3
"""Generate clean, symmetric TikZ subpanels for all 50 tight 4-loop cubic Feynman graphs."""
import json, math
from pathlib import Path

SRC = Path("/Users/libotao/.gemini/antigravity/scratch/holomorphic_5d_cs_paper/github_repo/deepen/o3/higher_orders/four-loop-graphs.json")
OUT = Path("/Users/libotao/.gemini/antigravity/scratch/holomorphic_5d_cs_paper/atlas_50_graphs.tex")

graphs = json.load(SRC.open())
assert len(graphs) == 50

def core_coords(idx):
    """Return {v: (x, y)} in [-1.05, 1.05]^2 for core vertices 0..5 by core family."""
    if idx <= 2:
        # SP Core (0..2): 0 left, 3 right, (1,2) top bubble, (5,4) bottom bubble, [0,3] chord
        return {
            0: (-0.95, 0.0),
            3: ( 0.95, 0.0),
            1: (-0.42, 0.62),
            2: ( 0.42, 0.62),
            5: (-0.42,-0.62),
            4: ( 0.42,-0.62),
        }
    elif idx <= 10:
        # Doubled-Kite Core (3..10): 5 center-top/hub, {0,1,2} arc, (3,4) bottom doubled bubble
        # Edges: [0,1],[0,4],[0,5],[1,2],[1,5],[2,3],[2,5],[3,4],[3,4]
        return {
            5: ( 0.0,  0.18),
            1: ( 0.0,  0.88),
            0: (-0.88, 0.30),
            2: ( 0.88, 0.30),
            4: (-0.48,-0.68),
            3: ( 0.48,-0.68),
        }
    elif idx <= 38:
        # Triangular Prism Core (11..38): inner triangle {0,2,3}, outer triangle {1,5,4}, rungs (0,1),(2,5),(3,4)
        pts = {}
        r_in, r_out = 0.38, 0.92
        angles = [90, 210, 330]
        for v, ang in zip([0, 3, 2], angles):
            rad = math.radians(ang)
            pts[v] = (r_in * math.cos(rad), r_in * math.sin(rad) - 0.06)
        for v, ang in zip([1, 4, 5], angles):
            rad = math.radians(ang)
            pts[v] = (r_out * math.cos(rad), r_out * math.sin(rad) - 0.06)
        return pts
    else:
        # Non-planar K_{3,3} Core (39..49): regular hexagon 0-1-2-3-4-5 + diagonals (0,3),(1,4),(2,5)
        pts = {}
        r = 0.88
        for k in range(6):
            rad = math.radians(90 - 60 * k)
            pts[k] = (r * math.cos(rad), r * math.sin(rad))
        return pts

def quad_bezier(p0, p1, pc, t):
    u = 1.0 - t
    return (
        u*u*p0[0] + 2*u*t*pc[0] + t*t*p1[0],
        u*u*p0[1] + 2*u*t*pc[1] + t*t*p1[1],
    )

def render_graph_cell(idx, g):
    coords = core_coords(idx)
    chains = g["chains"]
    aut = g["aut"]
    if idx <= 2:
        fam, col = "SP", "PineGreen!80!black"
    elif idx <= 10:
        fam, col = "Kite", "RoyalBlue!85!black"
    elif idx <= 38:
        fam, col = "Prism", "NavyBlue!85!black"
    else:
        fam, col = "$K_{3,3}$", "Plum!85!black"

    # Group chains by unordered core endpoints to detect parallel bubbles
    by_pair = {}
    for ci, ch in enumerate(chains):
        pair = (min(ch[0], ch[-1]), max(ch[0], ch[-1]))
        by_pair.setdefault(pair, []).append((ci, ch))

    lines = []
    lines.append(r"\begin{tikzpicture}[scale=0.88, baseline={(0,0)}]")
    lines.append(rf"  \draw[rounded corners=2.5pt, draw={col}!45, fill={col}!3, line width=0.45pt] (-1.28,-1.18) rectangle (1.28,1.26);")
    lines.append(rf"  \node[anchor=north west, font=\tiny\bfseries, text={col}, inner sep=1.5pt] at (-1.22,1.21) {{\#{idx:02d} ({fam})}};")
    lines.append(rf"  \node[anchor=north east, font=\tiny, text=black!75, inner sep=1.5pt] at (1.22,1.21) {{$|\mathrm{{Aut}}|={aut}$}};")

    ext_markers = []
    for pair, ch_list in by_pair.items():
        multi = len(ch_list)
        for rank, (ci, ch) in enumerate(ch_list):
            u, v = ch[0], ch[-1]
            p0, p1 = coords[u], coords[v]
            dx, dy = p1[0] - p0[0], p1[1] - p0[1]
            dist = math.hypot(dx, dy)
            nx, ny = -dy / dist, dx / dist  # unit normal
            # Offset for parallel edges (or slight curve on K3,3 diagonals so center doesn't triple-collide)
            offset = 0.0
            if multi == 2:
                # Orient normal consistently with pair order
                pu, pv = coords[pair[0]], coords[pair[1]]
                pdx, pdy = pv[0] - pu[0], pv[1] - pu[1]
                pnx, pny = -pdy / dist, pdx / dist
                offset = 0.24 if rank == 0 else -0.24
                pc = (0.5*(p0[0]+p1[0]) + offset*pnx, 0.5*(p0[1]+p1[1]) + offset*pny)
            elif idx >= 39 and pair in ((0,3), (1,4), (2,5)):
                # Slight bow on two of the 3 main diagonals in K3,3 so center is visually distinct
                bow = { (0,3): 0.0, (1,4): 0.14, (2,5): -0.14 }[pair]
                pc = (0.5*(p0[0]+p1[0]) + bow*nx, 0.5*(p0[1]+p1[1]) + bow*ny)
                offset = bow
            else:
                pc = (0.5*(p0[0]+p1[0]), 0.5*(p0[1]+p1[1]))

            if abs(offset) > 1e-6:
                lines.append(rf"  \draw[line width=0.6pt, draw=black!80] ({p0[0]:.3f},{p0[1]:.3f}) .. controls ({pc[0]:.3f},{pc[1]:.3f}) .. ({p1[0]:.3f},{p1[1]:.3f});")
            else:
                lines.append(rf"  \draw[line width=0.6pt, draw=black!80] ({p0[0]:.3f},{p0[1]:.3f}) -- ({p1[0]:.3f},{p1[1]:.3f});")

            # Subdivision bivalent external vertices along ch[1:-1]
            m = len(ch) - 1  # number of segments
            for k in range(1, m):
                t = k / m
                if abs(offset) > 1e-6:
                    # Quadratic Bezier with control point shifted so midpoint hits 0.5*offset:
                    # B(0.5) = 0.25 p0 + 0.5 pc + 0.25 p1 = mid + 0.5*offset*n
                    bx, by = quad_bezier(p0, p1, pc, t)
                    # Tangent at t: 2(1-t)(pc-p0) + 2t(p1-pc)
                    tx = 2*(1-t)*(pc[0]-p0[0]) + 2*t*(p1[0]-pc[0])
                    ty = 2*(1-t)*(pc[1]-p0[1]) + 2*t*(p1[1]-pc[1])
                else:
                    bx = p0[0] + t * dx
                    by = p0[1] + t * dy
                    tx, ty = dx, dy
                tlen = math.hypot(tx, ty) or 1.0
                # Outward normal for external leg stub (point away from origin when possible)
                enx, eny = -ty / tlen, tx / tlen
                if enx * bx + eny * by < -0.02:
                    enx, eny = -enx, -eny
                sx, sy = bx + 0.22 * enx, by + 0.22 * eny
                # Clamp stub inside box
                sx = max(-1.20, min(1.20, sx))
                sy = max(-1.10, min(0.96, sy))
                ext_markers.append((bx, by, sx, sy))

    for bx, by, sx, sy in ext_markers:
        lines.append(rf"  \draw[line width=0.65pt, draw=BrickRed!85!black] ({bx:.3f},{by:.3f}) -- ({sx:.3f},{sy:.3f});")
        lines.append(rf"  \filldraw[fill=BrickRed!85!black, draw=white, line width=0.25pt] ({bx:.3f},{by:.3f}) circle (1.35pt);")

    for v, (vx, vy) in coords.items():
        lines.append(rf"  \filldraw[fill=black!85, draw=white, line width=0.25pt] ({vx:.3f},{vy:.3f}) circle (1.45pt);")

    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines)

def build_panel(start, stop, fig_label, caption):
    out = []
    out.append(r"\begin{figure}[H]")
    out.append(r"\centering")
    out.append(r"\setlength{\tabcolsep}{2.2pt}")
    out.append(r"\renewcommand{\arraystretch}{1.05}")
    out.append(r"\begin{tabular}{@{}ccccc@{}}")
    for row_start in range(start, stop, 5):
        cells = [render_graph_cell(i, graphs[i]) for i in range(row_start, row_start + 5)]
        out.append(" &\n".join(cells) + r" \\[2pt]")
    out.append(r"\end{tabular}")
    out.append(rf"\caption{{{caption}}}")
    out.append(rf"\label{{{fig_label}}}")
    out.append(r"\end{figure}")
    return "\n".join(out)

tex = (
    "% Auto-generated exact 50-graph 4-loop Feynman atlas from four-loop-graphs.json\n"
    + build_panel(
        0, 25, "fig:atlas-00-24",
        r"Complete bare four-loop Feynman diagram atlas on $\mathbb{R}\times\mathbb{C}^2$, Part~I: Graphs \texttt{\#00}--\texttt{\#24} (comprising all $3$ series--parallel topologies \texttt{\#00}--\texttt{\#02}, all $8$ doubled-kite topologies \texttt{\#03}--\texttt{\#10}, and triangular-prism topologies \texttt{\#11}--\texttt{\#24}). Black nodes denote internal trivalent $f_{ABC}$ vertices ($6$ per graph); crimson nodes with outward stubs denote the $4$ bivalent external superfield insertions $\mathcal{A}$."
    )
    + "\n\n"
    + build_panel(
        25, 50, "fig:atlas-25-49",
        r"Complete bare four-loop Feynman diagram atlas on $\mathbb{R}\times\mathbb{C}^2$, Part~II: Graphs \texttt{\#25}--\texttt{\#49} (comprising the remaining $14$ triangular-prism topologies \texttt{\#25}--\texttt{\#38} and all $11$ non-planar bipartite $K_{3,3}$ topologies \texttt{\#39}--\texttt{\#49}, where antipodal chords are bowed slightly around the center to keep every crossing legible)."
    )
    + "\n"
)
OUT.write_text(tex)
print(f"Wrote {OUT} ({len(tex)} bytes)")
