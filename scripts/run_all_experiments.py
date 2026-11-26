"""Generate all 16 canonical HIRSCH-DIAM-CERT-1.0 certificates and 4 experiment JSONs.

Uses exact integer/Fraction arithmetic for all certificate generation and sweeps:
  - Exp 1: 16-certificate atlas & warm-up timing (4 per mode)
  - Exp 2: Superpolynomial Theorem 1.1/1.2 diameter vs Hirsch (n - d) & Kalai-Kleitman bounds
  - Exp 3: Exact rational Stage-2 smoothing deviation |B_{ij} - 1/d| vs smoothing parameter t
  - Exp 4: Exact Klee-Walkup wedge vertex/facet progression from pentagon/hexagon/3D prism to complementary 2d-facet polytopes
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import sys
import time
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.hirsch_certifier import PROTOCOL_ID, verify_certificate_dict
from src.hirsch_exact import (
    enumerate_simple_polytope_vertices,
    frac_div,
    klee_walkup_wedge_system,
    solve_linear_exact,
    stage1_doubly_stochastic,
    stage2_symmetric_smooth_packing,
)

CERT_DIR = ROOT / "certificates"
RES_DIR = ROOT / "results"


def _fp(f: Fraction) -> List[int]:
    return [f.numerator, f.denominator]


# ---------------------------------------------------------------------------
# Certificate Builders (4 per Mode = 16 Total)
# ---------------------------------------------------------------------------

def build_mode1_cert(
    cert_id: str,
    corners: List[Fraction],
    N: int,
    K: int,
    R: int,
    g: Fraction,
    d: int,
    E: int,
    d0: int,
    f0: int,
) -> Dict[str, Any]:
    # Build an exact rational zigzag path hitting z=0 at every corner a_k with left slope <= -1, right slope >= +1
    pts: List[List[List[int]]] = []
    # Start left of first corner with z = 2 * step
    step0 = Fraction(1, 2)
    pts.append([_fp(corners[0] - step0), _fp(Fraction(1))])  # slope to corner 0 is -2 <= -1
    for i, a in enumerate(corners):
        pts.append([_fp(a), _fp(Fraction(0))])
        if i + 1 < len(corners):
            mid = frac_div(a + corners[i + 1], Fraction(2))
            half = frac_div(corners[i + 1] - a, Fraction(2))
            # peak height = 2 * half => right slope from a is +2 >= 1, left slope into a_{i+1} is -2 <= -1!
            pts.append([_fp(mid), _fp(Fraction(2) * half)])
    pts.append([_fp(corners[-1] + step0), _fp(Fraction(1))])  # right slope +2 >= 1
    # Every interior waypoint (2*len(corners) - 1 of them) is a bend!
    bends = len(pts) - 2
    full_lb = R * (N ** K)
    gap_lb = Fraction(R * (N ** (K - 1))) * max(Fraction(0), Fraction(N) * g - Fraction(8))
    return {
        "protocol": PROTOCOL_ID,
        "certificate_id": cert_id,
        "mode": "convex_corner_bend_and_crossing_cert",
        "description": f"Secs. 1.2/2.1/2.2 {len(corners)} V-corners ({bends} breaks), Lem. 2.2 N={N}, K={K}, R={R}, Thm. 1.2 ({d}, {d+E})",
        "payload": {
            "corner_xs": [_fp(a) for a in corners],
            "path_points": pts,
            "claimed_breaks": bends,
            "N_param": N,
            "K_param": K,
            "R_param": R,
            "gap_g": _fp(g),
            "claimed_full_dist_lb": full_lb,
            "claimed_gap_dist_lb": _fp(gap_lb),
            "prod_d": d,
            "prod_E": E,
            "seed_d0": d0,
            "seed_f": f0,
        },
    }


def build_mode2_cert(
    cert_id: str,
    p: int,
    normals: List[List[Fraction]],
    rhs: List[Fraction],
    tip_v_idx: int,
    tip_w_idx: int,
    wedge_seq: List[int],
) -> Dict[str, Any]:
    cur_n = [list(r) for r in normals]
    cur_b = list(rhs)
    cur_d = p
    for f_idx in wedge_seq:
        cur_n, cur_b = klee_walkup_wedge_system(cur_n, cur_b, f_idx)
        cur_d += 1
    final_v, _ = enumerate_simple_polytope_vertices(cur_n, cur_b, cur_d)
    return {
        "protocol": PROTOCOL_ID,
        "certificate_id": cert_id,
        "mode": "klee_walkup_wedge_and_complement_cert",
        "description": f"Cor. 14.1 Klee-Walkup wedge ({p}D/{len(normals)}f -> {cur_d}D/{len(cur_n)}f, {len(final_v)} simple vertices)",
        "payload": {
            "base_dim": p,
            "base_normals": [[_fp(x) for x in r] for r in normals],
            "base_rhs": [_fp(b) for b in rhs],
            "tip_v_idx": tip_v_idx,
            "tip_w_idx": tip_w_idx,
            "wedge_facet_sequence": wedge_seq,
            "claimed_final_dim": cur_d,
            "claimed_final_facets": len(cur_n),
            "claimed_final_vertex_count": len(final_v),
        },
    }


def build_mode3_cert(
    cert_id: str,
    dm1: int,
    U: List[List[Fraction]],
    V: List[List[Fraction]],
    facets: List[Tuple[List[Fraction], Fraction, Fraction]],
    Ds: int,
) -> Dict[str, Any]:
    mid_pts = {tuple(frac_div(u[j] + v[j], Fraction(2)) for j in range(dm1)) for u in U for v in V}
    return {
        "protocol": PROTOCOL_ID,
        "certificate_id": cert_id,
        "mode": "schlegel_projective_and_prismatoid_cert",
        "description": f"Sec. 14.1 Cor. 14.2 {dm1+1}D prismatoid ({len(facets)} lateral facets, {len(mid_pts)} (U+V)/2 sums, Ds={Ds})",
        "payload": {
            "base_dim_d_minus_1": dm1,
            "T0_vertices": [[_fp(c) for c in pt] for pt in U],
            "T1_vertices": [[_fp(c) for c in pt] for pt in V],
            "nonbase_facets": [{"a": [_fp(x) for x in a], "b": _fp(b), "c": _fp(c)} for a, b, c in facets],
            "claimed_minkowski_vertex_count": len(mid_pts),
            "prismatoid_width_Ds": Ds,
            "claimed_section_dist": Ds - 2,
            "claimed_schlegel_depth": Ds - 1,
        },
    }


def build_mode4_cert(
    cert_id: str,
    A: List[List[Fraction]],
    M: Fraction,
    test_x: List[Fraction],
    odd_pad: bool,
    t_sm: Fraction,
    eps: Fraction,
) -> Dict[str, Any]:
    k = len(A)
    _, _, B0 = stage1_doubly_stochastic(A, M)
    d = 2 * k + (1 if odd_pad else 0)
    B = stage2_symmetric_smooth_packing(B0, odd_pad, t_sm)
    target = Fraction(1, d)
    max_dev = max(abs(B[i][j] - target) for i in range(d) for j in range(d))
    return {
        "protocol": PROTOCOL_ID,
        "certificate_id": cert_id,
        "mode": "symmetric_doubly_stochastic_packing_cert",
        "description": f"Cor. 14.3 two-stage packing ({k}x{k} -> symmetric {d}x{d} doubly stochastic, max dev {max_dev} < {eps})",
        "payload": {
            "k_dim": k,
            "A_matrix": [[_fp(v) for v in row] for row in A],
            "M_param": _fp(M),
            "test_x": [_fp(v) for v in test_x],
            "odd_pad": odd_pad,
            "t_smooth": _fp(t_sm),
            "epsilon_bound": _fp(eps),
            "claimed_d": d,
            "claimed_max_dev": _fp(max_dev),
        },
    }


def build_simplicial_antiprism(
    dm1: int,
    scale_v: Fraction = Fraction(-1, 3),
) -> Tuple[List[List[Fraction]], List[List[Fraction]], List[Tuple[List[Fraction], Fraction, Fraction]]]:
    """Build exact d-dimensional simplicial antiprism conv((U, 0) U (scale_v * U, 1))
    where U = {e_1, ..., e_{d-1}, -sum e_j} has barycenter 0 in int(conv(U)).
    Has 2^d - 2 simplicial lateral facets, every facet satisfies Schlegel strictness b + c > 0
    for scale_v = -1/3 and d in {3, 4, 5}, and exact dual-graph BFS width D_1 == d!
    """
    d = dm1 + 1
    U: List[List[Fraction]] = []
    for i in range(dm1):
        row = [Fraction(0)] * dm1
        row[i] = Fraction(1)
        U.append(row)
    U.append([Fraction(-1)] * dm1)
    V: List[List[Fraction]] = [[scale_v * c for c in u] for u in U]

    facets: List[Tuple[List[Fraction], Fraction, Fraction]] = []
    # Every nonempty proper subset I of {0..d-1} gives a simplicial facet spanning U[I] and V[J]
    for r in range(1, d):
        for I_tup in combinations(range(d), r):
            I_set = set(I_tup)
            # Linear functional a . u_i = 1 for i in I, a . u_j = -r / (d - r) for j in J
            # Since u_0..u_{dm1-1} are standard basis vectors, a[k] = target[k]!
            neg_val = Fraction(-r, d - r)
            a_vec = [Fraction(1) if k in I_set else neg_val for k in range(dm1)]
            # On V[j] (j not in I): a . (scale_v * u_j) + b = 1 => b = 1 - scale_v * (-r / (d - r))
            b_val = Fraction(1) - scale_v * neg_val
            c_val = Fraction(1)
            facets.append((a_vec, b_val, c_val))
    return U, V, facets


def main() -> None:
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    for old_cert in CERT_DIR.glob("*.cert.json"):
        old_cert.unlink()
    RES_DIR = ROOT / "results"
    RES_DIR.mkdir(parents=True, exist_ok=True)

    # --- Base polytopes for Mode 2 (Klee-Walkup wedges) ---
    # 1. Pentagon (2D, 5 facets -> 1 wedge -> 3D simple 6-facet complementary polytope!)
    pent_n = [
        [Fraction(-1), Fraction(0)],   # x >= 0
        [Fraction(0), Fraction(-1)],   # y >= 0
        [Fraction(1), Fraction(0)],    # x <= 2
        [Fraction(0), Fraction(1)],    # y <= 2
        [Fraction(1), Fraction(1)],    # x + y <= 3
    ]
    pent_b = [Fraction(0), Fraction(0), Fraction(2), Fraction(2), Fraction(3)]
    pent_v, _ = enumerate_simple_polytope_vertices(pent_n, pent_b, 2)
    # Tips: (0, 0) [facets 0,1] and (2, 1) [facets 2,4] -> uncovered facet is 3 (y <= 2)!
    p_iv = pent_v.index((Fraction(0), Fraction(0)))
    p_iw = pent_v.index((Fraction(2), Fraction(1)))

    # 2. Hexagon (2D, 6 facets -> 2 wedges -> 4D simple 8-facet complementary polytope!)
    hex_n = [
        [Fraction(-1), Fraction(0)],   # 0: x >= 0
        [Fraction(0), Fraction(-1)],   # 1: y >= 0
        [Fraction(1), Fraction(0)],    # 2: x <= 3
        [Fraction(0), Fraction(1)],    # 3: y <= 3
        [Fraction(1), Fraction(1)],    # 4: x + y <= 5
        [Fraction(-1), Fraction(1)],   # 5: -x + y <= 2
    ]
    hex_b = [Fraction(0), Fraction(0), Fraction(3), Fraction(3), Fraction(5), Fraction(2)]
    hex_v, _ = enumerate_simple_polytope_vertices(hex_n, hex_b, 2)
    # Tips: (0,0) [facets 0,1] and (3,2) [facets 2,4] -> uncovered facets are 3 and 5!
    # Note: after wedging facet 3 (split into 3,4), old facet 5 shifts to index 6!
    h_iv = hex_v.index((Fraction(0), Fraction(0)))
    h_iw = hex_v.index((Fraction(3), Fraction(2)))

    # 3. Heptagon (2D, 7 facets -> 3 wedges -> 5D simple 10-facet complementary polytope!)
    hep_n = hex_n + [[Fraction(1), Fraction(-1)]]  # 6: x - y <= 2
    hep_b = hex_b + [Fraction(2)]
    hep_v, _ = enumerate_simple_polytope_vertices(hep_n, hep_b, 2)
    # Tips: (0,0) [facets 0,1] and (3,2) [facets 2,4] -> uncovered facets: 3, 5, 6!
    # After wedging 3 -> 3,4; old 5->6, old 6->7. After wedging 6 -> 6,7; old 7->8!
    hp_iv = hep_v.index((Fraction(0), Fraction(0)))
    hp_iw = hep_v.index((Fraction(3), Fraction(2)))

    # 4. 3D Pentagonal Prism (3D, 7 facets -> 1 wedge -> 4D simple 8-facet complementary polytope!)
    prism_n = [r + [Fraction(0)] for r in pent_n] + [
        [Fraction(0), Fraction(0), Fraction(-1)],  # 5: z >= 0
        [Fraction(0), Fraction(0), Fraction(1)],   # 6: z <= 1
    ]
    prism_b = pent_b + [Fraction(0), Fraction(1)]
    prism_v, _ = enumerate_simple_polytope_vertices(prism_n, prism_b, 3)
    # Tips: (0,0,0) [facets 0,1,5] and (2,1,1) [facets 2,4,6] -> uncovered facet is 3 (y <= 2)!
    pr_iv = prism_v.index((Fraction(0), Fraction(0), Fraction(0)))
    pr_iw = prism_v.index((Fraction(2), Fraction(1), Fraction(1)))

    # --- Mode 3: Exact Prismatoids (U at z=0, V at z=1) ---
    # 3D triangular prismatoid (antiprism / prism): U = {(0,0),(2,0),(0,2)}, V = {(2,2),(0,2),(2,0)}?
    # Simpler: triangular prism U = V = {(0,0),(1,0),(0,1)} (3 lateral facets: -x<=0, -y<=0, x+y<=1; b=0, c>=0!)
    # Wait: viewpoint strictness needs b + c > 0! Shift U, V so (0,0) is strictly interior!
    # Take U = {(-1,-1), (2,-1), (-1,2)}, V = U (prism) or V = rotated {(1,1),(-2,1),(1,-2)} (octahedron-antiprism!).
    # Let's check triangular prism with U = V = {(-1,-1), (2,-1), (-1,2)}:
    #   Facet 1: -x <= 1 (a=[-1,0], b=0, c=1 -> b+c=1>0! Tight on 4 verts!)
    #   Facet 2: -y <= 1 (a=[0,-1], b=0, c=1 -> b+c=1>0! Tight on 4 verts!)
    #   Facet 3: x + y <= 1 (a=[1,1], b=0, c=1 -> b+c=1>0! Tight on 4 verts!)
    U3 = [[Fraction(-1), Fraction(-1)], [Fraction(2), Fraction(-1)], [Fraction(-1), Fraction(2)]]
    V3_prism = U3
    fac3_prism = [
        ([Fraction(-1), Fraction(0)], Fraction(0), Fraction(1)),
        ([Fraction(0), Fraction(-1)], Fraction(0), Fraction(1)),
        ([Fraction(1), Fraction(1)], Fraction(0), Fraction(1)),
    ]
    # Triangular antiprism (Figure 7 of preprint!): U = {(-1,-1),(2,-1),(-1,2)}, V = {(1,1),(-2,1),(1,-2)}
    # Check lateral facets of conv((U,0) U (V,1)):
    #   -x - y + 2z <= 2? At (-1,-1,0): 2<=2 (tight); at (2,-1,0): -1; at (-1,2,0): -1;
    #                     at (1,1,1): 0; at (-2,1,1): 3>2 (no).
    #   Exact antiprism: U = {(1,0), (-1,1), (0,-1)}? Better: U = {(2,0),(-1,2),(-1,-2)}?
    #   Even cleaner! U = {(1,0), (0,1), (-1,-1)}, V = {(-1,0), (0,-1), (1,1)}!
    #   Check face spanning u1=(1,0,0), u2=(0,1,0), v3=(1,1,1): plane x + y - z = 1!
    #     At U: (1,0,0)->1, (0,1,0)->1, (-1,-1,0)->-2 <= 1!
    #     At V: (-1,0,1)->-2, (0,-1,1)->-2, (1,1,1)->1 <= 1!
    #     b = -1, c = 1 -> b+c = 0! Wait: viewpoint (0,-1) has z=-1, so x+y-z = 1 at (0,-1)!
    #   Make V slightly scaled by 1/2: V = {(-1/2,0), (0,-1/2), (1/2,1/2)} (strictly inside cone from (0,-1)!):
    #     Plane through (1,0,0), (0,1,0), (1/2,1/2,1): x + y + 0*z = 1 -> b=0, c=1, b+c=1 > 0!!
    #     Check all 6 facets of triangular antiprism with U = {(1,0),(0,1),(-1,-1)}, V = {(-1/2,0),(0,-1/2),(1/2,1/2)}:
    #     1) x + y <= 1: U -> 1,1,-2; V -> -1/2,-1/2,1. Tight on 3 verts (u0,u1,v2)! b+c = 1 > 0!
    #     2) -2x + y + (1/2)z <= 1: U -> -2,1,1; V -> 1+1/2=3/2 > 1 -> adjust b!
    #     Solve exact planes for U={(1,0),(0,1),(-1,-1)}, V={(-1/2,0),(0,-1/2),(1/2,1/2)}:
    #       (u0,u1,v2): x + y <= 1 (b=0, c=1, b+c=1>0)
    #       (u1,u2,v0): -2x + y - (1/2)z <= 1? Wait: at u1=(0,1,0)->1; u2=(-1,-1,0)->1; v0=(-1/2,0,1)->1 - 1/2 = 1/2 -> b=0 gives -2(-1/2)+0=1!!
    #       OH!! Because V = -1/2 * U!! Every edge of U has exact normal with -1/2*U on the SAME vertical line:
    #       1) x + y <= 1 (u0,u1,v2; b=0, c=1)
    #       2) -2x + y <= 1 (u1,u2,v0; b=0, c=1)
    #       3) x - 2y <= 1 (u2,u0,v1; b=0, c=1)
    #       And the 3 upper-base triangles (v0,v1,u2), (v1,v2,u0), (v2,v0,u1):
    #       4) -x - y + (3/2)z <= 2? Check (u2=(-1,-1,0)->2, v0=(-1/2,0,1)->1/2+3/2=2, v1=(0,-1/2,1)->2; u0->-1, u1->-1, v2->-1+3/2=1/2 <= 2)! b+c = 7/2 > 0!!
    #       5) 2x - y + (3/2)z <= 2: (u0=(1,0,0)->2, v1=(0,-1/2,1)->1/2+3/2=2, v2=(1/2,1/2,1)->1/2+3/2=2; u1->-1, u2->-1, v0->-1+3/2=1/2 <= 2)! b+c = 7/2 > 0!!
    #       6) -x + 2y + (3/2)z <= 2: (u1=(0,1,0)->2, v2=(1/2,1/2,1)->2, v0=(-1/2,0,1)->2; others <= 2)! b+c = 7/2 > 0!!
    #     ALL 6 LATERAL FACETS OF THE TRIANGULAR ANTIPRISM (FIG. 7) ARE EXACT RATIONAL SIMPLICES WITH b+c > 0!!
    U3_anti = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)], [Fraction(-1), Fraction(-1)]]
    V3_anti = [[Fraction(-1, 2), Fraction(0)], [Fraction(0), Fraction(-1, 2)], [Fraction(1, 2), Fraction(1, 2)]]
    fac3_anti = [
        ([Fraction(1), Fraction(1)], Fraction(0), Fraction(1)),
        ([Fraction(-2), Fraction(1)], Fraction(0), Fraction(1)),
        ([Fraction(1), Fraction(-2)], Fraction(0), Fraction(1)),
        ([Fraction(-1), Fraction(-1)], Fraction(3, 2), Fraction(2)),
        ([Fraction(2), Fraction(-1)], Fraction(3, 2), Fraction(2)),
        ([Fraction(-1), Fraction(2)], Fraction(3, 2), Fraction(2)),
    ]

    # 4D and 5D Exact Simplicial Antiprisms (exact dual widths D_1 = 4 and D_1 = 5!)
    U4_anti, V4_anti, fac4_anti = build_simplicial_antiprism(3, Fraction(-1, 3))
    U5_anti, V5_anti, fac5_anti = build_simplicial_antiprism(4, Fraction(-1, 3))

    # --- Mode 4: Invertible Nonnegative Packing Matrices A (A >= 0, A*1 = 1) ---
    A2 = [[Fraction(3, 4), Fraction(1, 4)], [Fraction(1, 5), Fraction(4, 5)]]
    A3 = [
        [Fraction(1, 2), Fraction(1, 4), Fraction(1, 4)],
        [Fraction(1, 6), Fraction(2, 3), Fraction(1, 6)],
        [Fraction(1, 5), Fraction(1, 5), Fraction(3, 5)],
    ]
    A4 = [
        [Fraction(1, 2), Fraction(1, 6), Fraction(1, 6), Fraction(1, 6)],
        [Fraction(1, 8), Fraction(1, 2), Fraction(1, 4), Fraction(1, 8)],
        [Fraction(1, 10), Fraction(1, 5), Fraction(1, 2), Fraction(1, 5)],
        [Fraction(1, 4), Fraction(1, 8), Fraction(1, 8), Fraction(1, 2)],
    ]

    certs: List[Dict[str, Any]] = [
        # Mode 1 (4 certs)
        build_mode1_cert("m1_corners_3_n16_k2", [Fraction(1), Fraction(3), Fraction(5)], 16, 2, 1, Fraction(3, 4), 20, 20, 4, 8),
        build_mode1_cert("m1_corners_5_n32_k3", [Fraction(1), Fraction(2), Fraction(4), Fraction(7), Fraction(9)], 32, 3, 2, Fraction(1, 2), 40, 60, 8, 20),
        build_mode1_cert("m1_corners_7_n64_k4", [Fraction(i) for i in range(1, 8)], 64, 4, 3, Fraction(1, 4), 100, 200, 12, 40),
        build_mode1_cert("m1_corners_9_n128_k5", [Fraction(i, 2) for i in range(1, 10)], 128, 5, 4, Fraction(1, 8), 256, 512, 16, 80),
        # Mode 2 (4 certs)
        build_mode2_cert("m2_wedge_pentagon_2d_to_3d", 2, pent_n, pent_b, p_iv, p_iw, [3]),
        build_mode2_cert("m2_wedge_hexagon_2d_to_4d", 2, hex_n, hex_b, h_iv, h_iw, [3, 6]),
        build_mode2_cert("m2_wedge_heptagon_2d_to_5d", 2, hep_n, hep_b, hp_iv, hp_iw, [3, 6, 8]),
        build_mode2_cert("m2_wedge_prism_3d_to_4d", 3, prism_n, prism_b, pr_iv, pr_iw, [3]),
        # Mode 3 (4 certs) - exact BFS dual widths D_1 = 3, 2, 4, 5!
        build_mode3_cert("m3_schlegel_antiprism_3d_fig7", 2, U3_anti, V3_anti, fac3_anti, 3),
        build_mode3_cert("m3_schlegel_prism_3d", 2, U3, V3_prism, fac3_prism, 2),
        build_mode3_cert("m3_schlegel_antiprism_4d", 3, U4_anti, V4_anti, fac4_anti, 4),
        build_mode3_cert("m3_schlegel_antiprism_5d", 4, U5_anti, V5_anti, fac5_anti, 5),
        # Mode 4 (4 certs)
        build_mode4_cert("m4_packing_k2_even_d4", A2, Fraction(2), [Fraction(1, 2), Fraction(1, 3)], False, Fraction(5), Fraction(1, 20)),
        build_mode4_cert("m4_packing_k2_odd_d5", A2, Fraction(3), [Fraction(1, 4), Fraction(2, 3)], True, Fraction(10), Fraction(1, 50)),
        build_mode4_cert("m4_packing_k3_even_d6", A3, Fraction(2), [Fraction(1, 3), Fraction(1, 2), Fraction(1, 4)], False, Fraction(20), Fraction(1, 100)),
        build_mode4_cert("m4_packing_k4_odd_d9", A4, Fraction(3), [Fraction(1, 2), Fraction(1, 4), Fraction(1, 3), Fraction(1, 5)], True, Fraction(50), Fraction(1, 400)),
    ]

    atlas_rows: List[Dict[str, Any]] = []
    for cert in certs:
        cid = cert["certificate_id"]
        path = CERT_DIR / f"{cid}.cert.json"
        path.write_text(json.dumps(cert, indent=2) + "\n", encoding="utf-8", newline="\n")
        summary = verify_certificate_dict(cert)
        times_ms: List[float] = []
        for _ in range(7):
            t0 = time.perf_counter()
            verify_certificate_dict(cert)
            times_ms.append((time.perf_counter() - t0) * 1000.0)
        mean_ms = sum(times_ms) / len(times_ms)
        atlas_rows.append(
            {
                "certificate_id": cid,
                "mode": cert["mode"],
                "mean_ms": round(mean_ms, 3),
                "summary": summary,
            }
        )
        print(f"[OK] {cid:38s} | {mean_ms:7.3f} ms | {summary}")

    (RES_DIR / "exp1_certificate_atlas.json").write_text(
        json.dumps({"protocol": PROTOCOL_ID, "count": len(atlas_rows), "certificates": atlas_rows}, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )

    # Exp 2: Exact certified Mode 1 Lemma 2.2 full vs gap distance bounds across C01..C04
    exp2_rows = []
    for idx, cert in enumerate(certs[:4]):
        pl = cert["payload"]
        gap_f = Fraction(pl["claimed_gap_dist_lb"][0], pl["claimed_gap_dist_lb"][1])
        exp2_rows.append(
            {
                "cert_id": f"C{idx+1:02d}",
                "N": pl["N_param"],
                "K": pl["K_param"],
                "R": pl["R_param"],
                "full_dist_lb": pl["claimed_full_dist_lb"],
                "gap_dist_lb": float(gap_f),
                "padded_facets": pl["prod_d"] + pl["prod_E"],
            }
        )
    (RES_DIR / "exp2_crossing_growth.json").write_text(
        json.dumps({"rows": exp2_rows}, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    # Exp 3: Exact rational Stage-2 smoothing deviation |B_{ij} - 1/d| vs t in [1/4 .. 256] (d=6 from A3)
    _, _, B0_3 = stage1_doubly_stochastic(A3, Fraction(2))
    exp3_rows = []
    for t_val in [Fraction(1, 4), Fraction(1, 2), Fraction(1), Fraction(2), Fraction(4), Fraction(8), Fraction(16), Fraction(32), Fraction(64), Fraction(128)]:
        B = stage2_symmetric_smooth_packing(B0_3, False, t_val)
        max_dev = max(abs(B[i][j] - Fraction(1, 6)) for i in range(6) for j in range(6))
        exp3_rows.append(
            {
                "t_float": float(t_val),
                "max_dev_str": f"{max_dev.numerator}/{max_dev.denominator}",
                "max_dev_float": float(max_dev),
                "asymp_bound": float(Fraction(1, 6 * t_val)),
            }
        )
    (RES_DIR / "exp3_smoothing_decay.json").write_text(
        json.dumps({"rows": exp3_rows}, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    # Exp 4: Exact Klee-Walkup wedge progression (dim, facets, uncovered, vertices) for pentagon/hexagon/heptagon/prism
    exp4_rows = [
        {"name": "Pentagon (2D, 5f)", "wedges": 1, "final_d": 3, "final_f": 6, "final_v": atlas_rows[4]["summary"]["final_vertices"]},
        {"name": "Hexagon (2D, 6f)", "wedges": 2, "final_d": 4, "final_f": 8, "final_v": atlas_rows[5]["summary"]["final_vertices"]},
        {"name": "Heptagon (2D, 7f)", "wedges": 3, "final_d": 5, "final_f": 10, "final_v": atlas_rows[6]["summary"]["final_vertices"]},
        {"name": "Pent. Prism (3D, 7f)", "wedges": 1, "final_d": 4, "final_f": 8, "final_v": atlas_rows[7]["summary"]["final_vertices"]},
    ]
    (RES_DIR / "exp4_wedge_progression.json").write_text(
        json.dumps({"rows": exp4_rows}, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    print("All 16 certificates and 4 experiment JSON files generated successfully.")


if __name__ == "__main__":
    main()
