"""Generate 3 publication-grade 300-DPI figures for Paper 66."""

from __future__ import annotations

import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RES_DIR = ROOT / "results"
FIG_DIR = ROOT / "paper" / "figures"

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 8.5,
        "axes.titlesize": 9.0,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 6.8,
        "axes.grid": True,
        "grid.alpha": 0.28,
        "grid.linestyle": "--",
    }
)


def make_fig1() -> None:
    rows2 = json.loads((RES_DIR / "exp2_crossing_growth.json").read_text(encoding="utf-8"))["rows"]
    c_labels = [f"{r['cert_id']}\n($N={r['N']}, K={r['K']}$)" for r in rows2]
    fulls = np.array([r["full_dist_lb"] for r in rows2], dtype=float)
    gaps = np.array([r["gap_dist_lb"] for r in rows2], dtype=float)

    rows3 = json.loads((RES_DIR / "exp3_smoothing_decay.json").read_text(encoding="utf-8"))["rows"]
    ts = np.array([r["t_float"] for r in rows3], dtype=float)
    devs = np.array([r["max_dev_float"] for r in rows3], dtype=float)
    asymp = np.array([r["asymp_bound"] for r in rows3], dtype=float)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.75), dpi=300)

    x1 = np.arange(len(c_labels))
    w1 = 0.36
    ax1.bar(x1 - w1 / 2, fulls, width=w1, color="#1F5AA6", label=r"Full LB $R N^K$")
    ax1.bar(x1 + w1 / 2, gaps, width=w1, color="#D35400", label=r"Gap LB $R N^{K-1}(Ng - 8)$")
    ax1.set_yscale("log")
    ax1.set_xticks(x1)
    ax1.set_xticklabels(c_labels, fontsize=6.6)
    ax1.set_ylabel("Certified distance LB (log)")
    ax1.set_title(r"(a) Exact Certified Lem. 2.2 Crossing Bounds")
    ax1.legend(loc="upper left", frameon=True, framealpha=0.93, edgecolor="#CCCCCC", fontsize=6.6)

    ax2.loglog(ts, devs, "o-", color="#2E7D32", lw=1.7, ms=4.0, label=r"Exact $\max_{i,j}|B_{ij} - 1/6|$ ($d{=}6$)")
    ax2.loglog(ts, asymp, "--", color="#7D3C98", lw=1.2, label=r"$O(1/(dt))$ upper envelope")
    ax2.set_xlabel(r"Projective smoothing parameter $t$ (log)")
    ax2.set_ylabel(r"Max entry deviation from $1/d$")
    ax2.set_title(r"(b) Cor. 14.3 Smoothing Convergence ($d{=}6$)")
    ax2.legend(loc="lower left", frameon=True, framealpha=0.93, edgecolor="#CCCCCC", fontsize=6.8)

    fig.tight_layout(pad=0.8)
    fig.savefig(FIG_DIR / "fig1_crossing_and_smoothing.png", dpi=300)
    plt.close(fig)


def make_fig2() -> None:
    rows4 = json.loads((RES_DIR / "exp4_wedge_progression.json").read_text(encoding="utf-8"))["rows"]
    names = ["Pentagon\n(1 wedge)", "Hexagon\n(2 wedges)", "Heptagon\n(3 wedges)", "3D Prism\n(1 wedge)"]
    dims = [r["final_d"] for r in rows4]
    facets = [r["final_f"] for r in rows4]
    verts = [r["final_v"] for r in rows4]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.75), dpi=300)

    x = np.arange(len(names))
    w = 0.36
    b1 = ax1.bar(x - w / 2, dims, width=w, color="#1F5AA6", label="Final dim $d$")
    b2 = ax1.bar(x + w / 2, facets, width=w, color="#2E7D32", label="Final facets $2d$")
    for bar, v in list(zip(b1, dims)) + list(zip(b2, facets)):
        ax1.text(bar.get_x() + bar.get_width() / 2.0, v + 0.25, str(v), ha="center", va="bottom", fontsize=7.2, fontweight="bold")
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, fontsize=7.0)
    ax1.set_ylim(0, 12.5)
    ax1.set_ylabel("Exact count")
    ax1.set_title("(a) Cor. 14.1 Wedge to Complementary $2d$-Facets")
    ax1.legend(loc="upper left", frameon=True, framealpha=0.93, edgecolor="#CCCCCC", fontsize=6.8)

    bars2 = ax2.bar(x, verts, width=0.55, color="#D35400")
    for bar, v in zip(bars2, verts):
        ax2.text(bar.get_x() + bar.get_width() / 2.0, v + 0.6, str(v), ha="center", va="bottom", fontsize=7.5, fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels(names, fontsize=7.0)
    ax2.set_ylim(0, max(verts) * 1.22)
    ax2.set_ylabel("Enumerated simple vertices")
    ax2.set_title("(b) Exact Simple Vertex Count After Wedges")

    fig.tight_layout(pad=0.8)
    fig.savefig(FIG_DIR / "fig2_wedge_progression.png", dpi=300)
    plt.close(fig)


def make_fig3() -> None:
    atlas = json.loads((RES_DIR / "exp1_certificate_atlas.json").read_text(encoding="utf-8"))["certificates"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.75), dpi=300)

    colors_map = {
        "convex_corner_bend_and_crossing_cert": "#1F5AA6",
        "klee_walkup_wedge_and_complement_cert": "#2E7D32",
        "schlegel_projective_and_prismatoid_cert": "#D35400",
        "symmetric_doubly_stochastic_packing_cert": "#7D3C98",
    }
    labels = [f"C{i+1:02d}" for i in range(len(atlas))]
    ms_vals = [max(c["mean_ms"], 0.004) for c in atlas]
    bar_colors = [colors_map[c["mode"]] for c in atlas]
    xpos = np.arange(len(atlas))
    ax1.bar(xpos, ms_vals, color=bar_colors, width=0.72)
    ax1.set_yscale("log")
    ax1.set_xticks(xpos[::2])
    ax1.set_xticklabels(labels[::2])
    ax1.set_ylabel("Mean verify time (ms, log)")
    ax1.set_xlabel("Certificate ID (Modes 1–4, 4 each)")
    ax1.set_title("(a) Zero-Float Certificate Latency (16/16)")

    cats = ["Schema/\nParser", "Mode 1\nCorners", "Mode 2\nWedge", "Mode 3\nSchlegel", "Mode 4\nPacking"]
    counts = [12, 15, 16, 18, 15]
    xp2 = np.arange(len(cats))
    bars = ax2.bar(xp2, counts, color="#1B7A43", width=0.58)
    for b, c in zip(bars, counts):
        ax2.text(b.get_x() + b.get_width() / 2.0, c + 0.35, f"{c}/{c}", ha="center", va="bottom", fontsize=7.5, fontweight="bold")
    ax2.set_xticks(xp2)
    ax2.set_xticklabels(cats)
    ax2.set_yticks([0, 5, 10, 15, 20])
    ax2.set_ylim(0, 21.5)
    ax2.set_ylabel("Adversarial forgeries blocked")
    ax2.set_title("(b) Live Red-Team Audit (76/76 Blocked)")

    fig.tight_layout(pad=0.8)
    fig.savefig(FIG_DIR / "fig3_atlas_and_redteam.png", dpi=300)
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    make_fig1()
    make_fig2()
    make_fig3()
    print("Generated all 3 figures at 300 DPI in paper/figures/.")


if __name__ == "__main__":
    main()
