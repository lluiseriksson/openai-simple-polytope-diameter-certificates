"""Zero-float 4-mode certificate auditor for HIRSCH-DIAM-CERT-1.0.

Strictly parses JSON without floats (rejecting IEEE-754 decimals, NaN/Infinity,
and duplicate keys via object_pairs_hook) and verifies four canonical modes
aligned with OpenAI's November 24, 2026 preprint "Superpolynomial graph diameter
of simple polytopes":
  1. convex_corner_bend_and_crossing_cert      (Secs. 1.2, 2.1, 2.2, Lem. 2.2, Thm. 1.2)
  2. klee_walkup_wedge_and_complement_cert     (Sec. 14, Cor. 14.1 exact simple wedge)
  3. schlegel_projective_and_prismatoid_cert   (Sec. 14.1, Cor. 14.2 Schlegel & (U+V)/2)
  4. symmetric_doubly_stochastic_packing_cert  (Sec. 14.2, Cor. 14.3 two-stage packing)

NO floating-point operations, true division '/', or math imports appear anywhere in this module.
"""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Tuple

from src.hirsch_exact import (
    CertificateError,
    affine_rank_exact,
    compute_prismatoid_dual_width,
    enumerate_simple_polytope_vertices,
    frac_div,
    gcd_int,
    klee_walkup_wedge_system,
    mat_det_exact,
    origin_strictly_inside_simplex,
    stage1_doubly_stochastic,
    stage2_symmetric_smooth_packing,
    verify_simple_polytope_bounded_and_irredundant,
)

PROTOCOL_ID = "HIRSCH-DIAM-CERT-1.0"
VALID_MODES = frozenset(
    {
        "convex_corner_bend_and_crossing_cert",
        "klee_walkup_wedge_and_complement_cert",
        "schlegel_projective_and_prismatoid_cert",
        "symmetric_doubly_stochastic_packing_cert",
    }
)


def _reject_float(val: str) -> Any:
    raise CertificateError(f"Floating-point literal rejected by zero-float policy: {val!r}")


def _reject_constant(val: str) -> Any:
    raise CertificateError(f"Non-standard JSON constant rejected: {val!r}")


def _no_duplicate_keys(pairs: List[Tuple[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for k, v in pairs:
        if k in out:
            raise CertificateError(f"Duplicate JSON key rejected: {k!r}")
        out[k] = v
    return out


def load_strict_json_str(raw: str) -> Dict[str, Any]:
    try:
        data = json.loads(
            raw,
            parse_float=_reject_float,
            parse_constant=_reject_constant,
            object_pairs_hook=_no_duplicate_keys,
        )
    except CertificateError:
        raise
    except Exception as exc:
        raise CertificateError(f"Malformed JSON: {exc}") from exc
    if not isinstance(data, dict):
        raise CertificateError("Top-level certificate payload must be a JSON object")
    return data


def _req_int(obj: Mapping[str, Any], key: str) -> int:
    if key not in obj:
        raise CertificateError(f"Missing required integer field {key!r}")
    v = obj[key]
    if isinstance(v, bool) or not isinstance(v, int):
        raise CertificateError(f"Field {key!r} must be an exact integer")
    return v


def _req_str(obj: Mapping[str, Any], key: str) -> str:
    if key not in obj:
        raise CertificateError(f"Missing required string field {key!r}")
    v = obj[key]
    if not isinstance(v, str) or not v.strip():
        raise CertificateError(f"Field {key!r} must be a non-empty string")
    return v


def _req_int_list(obj: Mapping[str, Any], key: str, min_len: int = 1, max_len: int = 64) -> List[int]:
    if key not in obj:
        raise CertificateError(f"Missing required integer list {key!r}")
    v = obj[key]
    if not isinstance(v, list) or not (min_len <= len(v) <= max_len):
        raise CertificateError(f"Field {key!r} must be a list of length in [{min_len}, {max_len}]")
    out: List[int] = []
    for item in v:
        if isinstance(item, bool) or not isinstance(item, int):
            raise CertificateError(f"Field {key!r} contains non-integer entry")
        out.append(item)
    return out


def _parse_frac_pair(pair: Any, label: str) -> Fraction:
    if not isinstance(pair, list) or len(pair) != 2:
        raise CertificateError(f"{label} must be [num, den]")
    if any(isinstance(x, bool) or not isinstance(x, int) for x in pair) or pair[1] <= 0:
        raise CertificateError(f"{label} must have integer numerator and positive integer denominator")
    if gcd_int(pair[0], pair[1]) != 1:
        raise CertificateError(f"{label} rational {pair} must be in lowest terms (gcd == 1)")
    return Fraction(pair[0], pair[1])


def _req_frac(obj: Mapping[str, Any], key: str) -> Fraction:
    if key not in obj:
        raise CertificateError(f"Missing required rational field {key!r}")
    return _parse_frac_pair(obj[key], key)


def _check_keys(obj: Mapping[str, Any], expected: set[str], label: str) -> None:
    keys = set(obj.keys())
    missing = expected - keys
    extra = keys - expected
    if missing:
        raise CertificateError(f"Missing keys in {label}: {sorted(missing)}")
    if extra:
        raise CertificateError(f"Unexpected keys in {label}: {sorted(extra)}")


# ---------------------------------------------------------------------------
# Mode 1: Convex Corner Bends, Lemma 2.2 Distance & Thm. 1.2 Product Counts
# ---------------------------------------------------------------------------

def _verify_mode1(payload: Mapping[str, Any]) -> Dict[str, Any]:
    expected = {
        "corner_xs",
        "path_points",
        "claimed_breaks",
        "N_param",
        "K_param",
        "R_param",
        "gap_g",
        "claimed_full_dist_lb",
        "claimed_gap_dist_lb",
        "prod_d",
        "prod_E",
        "seed_d0",
        "seed_f",
    }
    _check_keys(payload, expected, "Mode 1 payload")

    # 1. Verify Figure 1 convex-corner lower bound z_i >= |x - a_i| forces an interior break at every corner
    raw_corners = payload.get("corner_xs")
    if not isinstance(raw_corners, list) or not (2 <= len(raw_corners) <= 12):
        raise CertificateError("corner_xs must contain 2..12 strictly increasing rationals")
    corners = [_parse_frac_pair(p, f"corner_xs[{i}]") for i, p in enumerate(raw_corners)]
    for i in range(len(corners) - 1):
        if corners[i] >= corners[i + 1]:
            raise CertificateError("corner_xs must be strictly increasing")

    raw_pts = payload.get("path_points")
    if not isinstance(raw_pts, list) or not (len(corners) + 2 <= len(raw_pts) <= 32):
        raise CertificateError("path_points must have at least len(corner_xs)+2 waypoints")
    pts: List[Tuple[Fraction, Fraction]] = []
    for i, pt in enumerate(raw_pts):
        if not isinstance(pt, list) or len(pt) != 2:
            raise CertificateError(f"path_points[{i}] must be [x_frac, z_frac]")
        pts.append((_parse_frac_pair(pt[0], f"pt[{i}].x"), _parse_frac_pair(pt[1], f"pt[{i}].z")))

    # Verify x strictly increases across waypoints, starts < corners[0], ends > corners[-1],
    # satisfies z >= max_k |x - a_k|? Wait: Section 1.2 uses independent coordinates z_k >= |x - a_k|,
    # or equivalently a piecewise convex envelope: at each corner a_k, the path hits z_k = 0!
    # Here we verify exact single-coordinate V-corner feasibility at each visited corner a_k:
    # if waypoint j has x_j == a_k and z_j == 0, then left slope (z_j - z_{j-1})/(x_j - x_{j-1}) <= -1
    # and right slope (z_{j+1} - z_j)/(x_{j+1} - x_j) >= 1, proving a mandatory bend (slope jump >= 2)!
    for i in range(len(pts) - 1):
        if pts[i][0] >= pts[i + 1][0]:
            raise CertificateError("path_points x-coordinates must be strictly increasing")
    if not (pts[0][0] < corners[0] and pts[-1][0] > corners[-1]):
        raise CertificateError("path_points must strictly span all corners")
    # Enforce piecewise V-corner lower envelope z >= min_k |x - a_k| (and thus z >= 0) at every waypoint
    for idx, (px, pz) in enumerate(pts):
        min_v = min(abs(px - a) for a in corners)
        if pz < min_v:
            raise CertificateError(f"path_points[{idx}] violates V-corner lower envelope z >= min_k |x - a_k|")

    corner_set = set(corners)
    bends = 0
    corners_hit = 0
    for i in range(1, len(pts) - 1):
        x_prev, z_prev = pts[i - 1]
        x_cur, z_cur = pts[i]
        x_next, z_next = pts[i + 1]
        s_left = frac_div(z_cur - z_prev, x_cur - x_prev)
        s_right = frac_div(z_next - z_cur, x_next - x_cur)
        if x_cur in corner_set:
            if z_cur != Fraction(0):
                raise CertificateError(f"Corner contact at x={x_cur} must attain z=0")
            if s_left > Fraction(-1) or s_right < Fraction(1):
                raise CertificateError(f"Corner {x_cur} violates V-bound z >= |x - a| (slopes {s_left}, {s_right})")
            corners_hit += 1
        if s_left != s_right:
            bends += 1

    if corners_hit != len(corners):
        raise CertificateError(f"Path hit {corners_hit} != {len(corners)} prescribed corners")
    claimed_breaks = _req_int(payload, "claimed_breaks")
    if claimed_breaks != bends or bends < len(corners):
        raise CertificateError(f"claimed_breaks={claimed_breaks} disagrees with exact {bends} (>= {len(corners)})")

    # 2. Lemma 2.2 exact distance bounds R * N^K and R * N^{K-1} * max(0, N*g - 8)
    N = _req_int(payload, "N_param")
    K = _req_int(payload, "K_param")
    R = _req_int(payload, "R_param")
    g = _req_frac(payload, "gap_g")
    if not (16 <= N <= 1024 and 2 <= K <= 6 and 1 <= R <= 1000):
        raise CertificateError("Lemma 2.2 parameters (N, K, R) out of bounded range")
    if not (Fraction(0) < g <= Fraction(1)):
        raise CertificateError("gap_g must lie in (0, 1]")

    full_lb = R * (N ** K)
    gap_factor = max(Fraction(0), Fraction(N) * g - Fraction(8))
    gap_lb = Fraction(R * (N ** (K - 1))) * gap_factor
    if _req_int(payload, "claimed_full_dist_lb") != full_lb:
        raise CertificateError(f"claimed_full_dist_lb disagrees with R*N^K={full_lb}")
    if _req_frac(payload, "claimed_gap_dist_lb") != gap_lb:
        raise CertificateError(f"claimed_gap_dist_lb disagrees with exact {gap_lb}")

    # 3. Theorem 1.2 product padding exact accounting (dimension d, facets d + E)
    d = _req_int(payload, "prod_d")
    E = _req_int(payload, "prod_E")
    d0 = _req_int(payload, "seed_d0")
    f0 = _req_int(payload, "seed_f")
    if not (d0 >= 2 and d0 + 2 <= d <= 100_000 and f0 >= d0 + 1 and E + d0 - f0 + 1 >= 3):
        raise CertificateError("Theorem 1.2 product padding dimensions/facets invalid")
    simplex_dim = d - d0 - 2
    simplex_facets = simplex_dim + 1
    poly_sides = E + d0 - f0 + 1
    if (d0 + simplex_dim + 2) != d or (f0 + simplex_facets + poly_sides) != (d + E):
        raise CertificateError("Theorem 1.2 product dimension/facet sum mismatch")

    return {
        "mode": "convex_corner_bend_and_crossing_cert",
        "corners_verified": corners_hit,
        "interior_breaks": bends,
        "full_dist_lb": full_lb,
        "gap_dist_lb": f"{gap_lb.numerator}/{gap_lb.denominator}",
        "padded_facets": d + E,
    }


# ---------------------------------------------------------------------------
# Mode 2: Exact Rational Klee-Walkup Wedge & Complementary Vertices (Cor. 14.1)
# ---------------------------------------------------------------------------

def _verify_mode2(payload: Mapping[str, Any]) -> Dict[str, Any]:
    expected = {
        "base_dim",
        "base_normals",
        "base_rhs",
        "tip_v_idx",
        "tip_w_idx",
        "wedge_facet_sequence",
        "claimed_final_dim",
        "claimed_final_facets",
        "claimed_final_vertex_count",
    }
    _check_keys(payload, expected, "Mode 2 payload")

    p = _req_int(payload, "base_dim")
    if not (2 <= p <= 3):
        raise CertificateError("base_dim must be 2 or 3 for exact exhaustive vertex-lattice wedge audit")

    raw_normals = payload.get("base_normals")
    raw_rhs = payload.get("base_rhs")
    if not isinstance(raw_normals, list) or not isinstance(raw_rhs, list) or len(raw_normals) != len(raw_rhs):
        raise CertificateError("base_normals and base_rhs must be lists of equal length")
    m = len(raw_normals)
    if not (2 * p + 1 <= m <= 2 * p + 3):
        raise CertificateError(f"base polytope must have m in [{2*p+1}, {2*p+3}] facets (uncovered facets 1..3)")

    normals: List[List[Fraction]] = []
    for i, row in enumerate(raw_normals):
        if not isinstance(row, list) or len(row) != p:
            raise CertificateError(f"base_normals[{i}] must have length {p}")
        parsed_row = [_parse_frac_pair(c, f"base_normals[{i}][{j}]") for j, c in enumerate(row)]
        if all(c == 0 for c in parsed_row):
            raise CertificateError(f"base_normals[{i}] cannot be the zero vector")
        normals.append(parsed_row)
    rhs = [_parse_frac_pair(b, f"base_rhs[{i}]") for i, b in enumerate(raw_rhs)]

    verts, tight_sets = enumerate_simple_polytope_vertices(normals, rhs, p)
    verify_simple_polytope_bounded_and_irredundant(verts, tight_sets, m, p)

    iv = _req_int(payload, "tip_v_idx")
    iw = _req_int(payload, "tip_w_idx")
    if not (0 <= iv < len(verts) and 0 <= iw < len(verts)) or iv == iw:
        raise CertificateError("tip_v_idx and tip_w_idx must be distinct valid base vertex indices")
    v_pt = list(verts[iv])
    w_pt = list(verts[iw])
    tv = set(tight_sets[iv])
    tw = set(tight_sets[iw])
    if tv & tw:
        raise CertificateError("Base tips v and w share a facet; cannot be disjoint tips")

    wedge_seq = _req_int_list(payload, "wedge_facet_sequence", min_len=1, max_len=3)
    if len(wedge_seq) != m - 2 * p:
        raise CertificateError(f"wedge_facet_sequence length {len(wedge_seq)} != uncovered facets {m - 2 * p}")

    cur_dim = p
    cur_normals = normals
    cur_rhs = rhs
    cur_v = v_pt
    cur_w = w_pt
    for step, f_idx in enumerate(wedge_seq):
        if not (0 <= f_idx < len(cur_normals)):
            raise CertificateError(f"wedge_facet_sequence[{step}]={f_idx} out of range")
        slack_v = cur_rhs[f_idx] - sum(cur_normals[f_idx][j] * cur_v[j] for j in range(cur_dim))
        slack_w = cur_rhs[f_idx] - sum(cur_normals[f_idx][j] * cur_w[j] for j in range(cur_dim))
        if slack_v <= 0 or slack_w <= 0:
            raise CertificateError(f"Wedge facet {f_idx} at step {step} is not strictly uncovered by both tips")
        cur_normals, cur_rhs = klee_walkup_wedge_system(cur_normals, cur_rhs, f_idx)
        # Lift v to bottom sheet (t=0) and w to top sheet (t=slack_w)
        cur_v = cur_v + [Fraction(0)]
        cur_w = cur_w + [slack_w]
        cur_dim += 1

    final_verts, final_tights = enumerate_simple_polytope_vertices(cur_normals, cur_rhs, cur_dim)
    verify_simple_polytope_bounded_and_irredundant(final_verts, final_tights, len(cur_normals), cur_dim)
    if len(cur_normals) != 2 * cur_dim:
        raise CertificateError("Final wedged polytope does not have exactly 2*d facets")
    if tuple(cur_v) not in final_verts or tuple(cur_w) not in final_verts:
        raise CertificateError("Lifted tips are not vertices of final wedged polytope")

    ft_v = set(final_tights[final_verts.index(tuple(cur_v))])
    ft_w = set(final_tights[final_verts.index(tuple(cur_w))])
    if (ft_v & ft_w) or (ft_v | ft_w) != set(range(2 * cur_dim)):
        raise CertificateError("Final lifted tips are NOT complementary across all 2*d facets")

    if _req_int(payload, "claimed_final_dim") != cur_dim:
        raise CertificateError("claimed_final_dim mismatch")
    if _req_int(payload, "claimed_final_facets") != len(cur_normals):
        raise CertificateError("claimed_final_facets mismatch")
    if _req_int(payload, "claimed_final_vertex_count") != len(final_verts):
        raise CertificateError(f"claimed_final_vertex_count != exact {len(final_verts)}")

    return {
        "mode": "klee_walkup_wedge_and_complement_cert",
        "base_dim": p,
        "wedges_applied": len(wedge_seq),
        "final_dim": cur_dim,
        "final_facets": len(cur_normals),
        "final_vertices": len(final_verts),
    }


# ---------------------------------------------------------------------------
# Mode 3: Schlegel Projective Map & Prismatoid Middle Section (Sec. 14.1)
# ---------------------------------------------------------------------------

def _verify_mode3(payload: Mapping[str, Any]) -> Dict[str, Any]:
    expected = {
        "base_dim_d_minus_1",
        "T0_vertices",
        "T1_vertices",
        "nonbase_facets",
        "claimed_minkowski_vertex_count",
        "prismatoid_width_Ds",
        "claimed_section_dist",
        "claimed_schlegel_depth",
    }
    _check_keys(payload, expected, "Mode 3 payload")

    dm1 = _req_int(payload, "base_dim_d_minus_1")
    if not (2 <= dm1 <= 4):
        raise CertificateError("base_dim_d_minus_1 must lie in [2, 4]")
    d = dm1 + 1

    def _parse_verts(key: str) -> List[Tuple[Fraction, ...]]:
        raw = payload.get(key)
        if not isinstance(raw, list) or len(raw) != d:
            raise CertificateError(f"{key} must contain exactly d={d} rational simplex vertices")
        out: List[Tuple[Fraction, ...]] = []
        for i, v in enumerate(raw):
            if not isinstance(v, list) or len(v) != dm1:
                raise CertificateError(f"{key}[{i}] must have length {dm1}")
            out.append(tuple(_parse_frac_pair(c, f"{key}[{i}][{j}]") for j, c in enumerate(v)))
        if len(set(out)) != d:
            raise CertificateError(f"{key} vertices must be distinct")
        # Check affine independence via determinant of [(v_i - v_0)]
        mat = [[out[i][j] - out[0][j] for j in range(dm1)] for i in range(1, d)]
        if mat_det_exact(mat) == 0:
            raise CertificateError(f"{key} does not span a nondegenerate {dm1}-simplex")
        return out

    U = _parse_verts("T0_vertices")
    V = _parse_verts("T1_vertices")
    if not origin_strictly_inside_simplex(U):
        raise CertificateError("Origin 0 must lie strictly in the interior of T0_vertices for Schlegel projection")

    # Verify middle section (U + V)/2 pairwise sums count
    mid_pts = {tuple(frac_div(u[j] + v[j], Fraction(2)) for j in range(dm1)) for u in U for v in V}
    if len(mid_pts) > d * d:
        raise CertificateError("Minkowski sum exceeded d^2 pairwise vertex sums")
    if _req_int(payload, "claimed_minkowski_vertex_count") != len(mid_pts):
        raise CertificateError(f"claimed_minkowski_vertex_count != exact {len(mid_pts)}")

    # Verify Schlegel projective transformation (x, z) -> (X, Z) = (x/(1+z), -z/(1+z))
    # on every non-base facet a . x + b * z <= c:
    #   (1) Non-zero horizontal normal a (excludes T0 and T1), canonical distinctness
    #   (2) Strictness at viewpoint (0, -1): a.0 + b*(-1) < c <=> b + c > 0
    #   (3) Feasibility on ALL base vertices (u, 0) and (v, 1) & exact affine rank d - 1 of tight vertices
    #   (4) Exact transformed inequality a . X - (b + c) * Z <= c holds on all mapped images,
    #       and every T1 vertex (v, 1) maps to uniform height Z = -1/2
    #   (5) Exact facet completeness and dual-graph BFS width D_1 == prismatoid_width_Ds!
    raw_facets = payload.get("nonbase_facets")
    if not isinstance(raw_facets, list) or not (3 <= len(raw_facets) <= 32):
        raise CertificateError("nonbase_facets must contain 3..32 non-base inequalities")

    all_verts_d: List[Tuple[Fraction, ...]] = [tuple(list(u) + [Fraction(0)]) for u in U] + [
        tuple(list(v) + [Fraction(1)]) for v in V
    ]
    seen_canon: set[Tuple[Fraction, ...]] = set()
    facet_tight_sets: List[Tuple[int, ...]] = [tuple(range(d)), tuple(range(d, 2 * d))]

    for idx, fac in enumerate(raw_facets):
        if not isinstance(fac, dict):
            raise CertificateError(f"nonbase_facets[{idx}] must be a JSON object")
        _check_keys(fac, {"a", "b", "c"}, f"nonbase_facets[{idx}]")
        raw_a = fac.get("a")
        if not isinstance(raw_a, list) or len(raw_a) != dm1:
            raise CertificateError(f"nonbase_facets[{idx}].a must have length {dm1}")
        a_vec = [_parse_frac_pair(val, f"fac[{idx}].a[{j}]") for j, val in enumerate(raw_a)]
        b_val = _parse_frac_pair(fac["b"], f"fac[{idx}].b")
        c_val = _parse_frac_pair(fac["c"], f"fac[{idx}].c")
        if all(x == 0 for x in a_vec):
            raise CertificateError(f"nonbase_facets[{idx}] has zero horizontal normal (base facet T0/T1 rejected)")
        # Canonicalize by first nonzero entry of a_vec (preserving sign of positive scaling via |first_nz|)
        first_nz = next(abs(x) for x in a_vec if x != 0)
        canon = tuple([frac_div(x, first_nz) for x in a_vec] + [frac_div(b_val, first_nz), frac_div(c_val, first_nz)])
        if canon in seen_canon:
            raise CertificateError(f"Duplicate non-base facet rejected at nonbase_facets[{idx}]")
        seen_canon.add(canon)

        if b_val + c_val <= 0:
            raise CertificateError(f"nonbase_facets[{idx}] violates Schlegel viewpoint strictness b + c > 0")
        tight_indices: List[int] = []
        for u_idx, u in enumerate(U):
            lhs = sum(a_vec[j] * u[j] for j in range(dm1))
            if lhs > c_val:
                raise CertificateError(f"nonbase_facets[{idx}] violated by T0 vertex {u}")
            if lhs == c_val:
                tight_indices.append(u_idx)
        for v_idx, v in enumerate(V):
            lhs = sum(a_vec[j] * v[j] for j in range(dm1)) + b_val
            if lhs > c_val:
                raise CertificateError(f"nonbase_facets[{idx}] violated by T1 vertex {v}")
            if lhs == c_val:
                tight_indices.append(d + v_idx)
            # Image of (v, 1) is (v/2, -1/2)
            X_img = [frac_div(v[j], Fraction(2)) for j in range(dm1)]
            Z_img = Fraction(-1, 2)
            lhs_img = sum(a_vec[j] * X_img[j] for j in range(dm1)) - (b_val + c_val) * Z_img
            # Exact algebraic identity: lhs_img - c_val == (lhs - c_val) / 2
            if lhs_img - c_val != frac_div(lhs - c_val, Fraction(2)):
                raise CertificateError("Schlegel projective slack identity failed on T1 vertex")
            if lhs_img > c_val:
                raise CertificateError("Schlegel image of T1 vertex violated transformed facet")
        tight_pts = [all_verts_d[ti] for ti in tight_indices]
        if len(tight_indices) < d or affine_rank_exact(tight_pts) != d - 1:
            raise CertificateError(
                f"nonbase_facets[{idx}] is not a supporting (d-1)-facet (tight={len(tight_indices)}, rank<{d-1})"
            )
        facet_tight_sets.append(tuple(sorted(tight_indices)))

    actual_Ds = compute_prismatoid_dual_width(all_verts_d, facet_tight_sets, d)
    Ds = _req_int(payload, "prismatoid_width_Ds")
    if Ds != actual_Ds:
        raise CertificateError(f"prismatoid_width_Ds={Ds} disagrees with exact dual BFS width {actual_Ds}")
    if _req_int(payload, "claimed_section_dist") != Ds - 2:
        raise CertificateError(f"claimed_section_dist != Ds - 2 = {Ds - 2}")
    if _req_int(payload, "claimed_schlegel_depth") != Ds - 1:
        raise CertificateError(f"claimed_schlegel_depth != Ds - 1 = {Ds - 1}")

    return {
        "mode": "schlegel_projective_and_prismatoid_cert",
        "prismatoid_dim": d,
        "facets_verified": len(raw_facets),
        "minkowski_pair_sums": len(mid_pts),
        "actual_width_Ds": actual_Ds,
        "section_dist": Ds - 2,
        "schlegel_depth": Ds - 1,
    }


# ---------------------------------------------------------------------------
# Mode 4: Corollary 14.3 Two-Stage Symmetric Doubly-Stochastic Packing
# ---------------------------------------------------------------------------

def _verify_mode4(payload: Mapping[str, Any]) -> Dict[str, Any]:
    expected = {
        "k_dim",
        "A_matrix",
        "M_param",
        "test_x",
        "odd_pad",
        "t_smooth",
        "epsilon_bound",
        "claimed_d",
        "claimed_max_dev",
    }
    _check_keys(payload, expected, "Mode 4 payload")

    k = _req_int(payload, "k_dim")
    if not (2 <= k <= 4):
        raise CertificateError("k_dim must lie in [2, 4]")
    odd_pad = payload.get("odd_pad")
    if not isinstance(odd_pad, bool):
        raise CertificateError("odd_pad must be a boolean")

    raw_A = payload.get("A_matrix")
    if not isinstance(raw_A, list) or len(raw_A) != k:
        raise CertificateError(f"A_matrix must be {k}x{k}")
    A: List[List[Fraction]] = []
    for i, row in enumerate(raw_A):
        if not isinstance(row, list) or len(row) != k:
            raise CertificateError(f"A_matrix[{i}] must have length {k}")
        A.append([_parse_frac_pair(v, f"A[{i}][{j}]") for j, v in enumerate(row)])

    # Check A >= 0 (packing LP), A * 1 = 1, and A invertible
    for i in range(k):
        if any(v < 0 for v in A[i]):
            raise CertificateError(f"A_matrix row {i} has negative entry (violates packing LP A >= 0)")
        if sum(A[i]) != Fraction(1):
            raise CertificateError(f"A_matrix row {i} sum != 1 (violates A*1 = 1)")
    if mat_det_exact(A) == 0:
        raise CertificateError("A_matrix must be invertible (det != 0)")

    M = _req_frac(payload, "M_param")
    if M <= 0:
        raise CertificateError("M_param must be positive")
    c_vec, h, B0 = stage1_doubly_stochastic(A, M)

    if any(cj <= 0 for cj in c_vec):
        raise CertificateError("Stage-1 vector c must be strictly positive (increase M_param)")
    if sum(c_vec) != h - Fraction(1):
        raise CertificateError("Stage-1 identity c^T 1 = h - 1 failed")
    for i in range(k):
        if sum(B0[i]) != Fraction(1) or sum(B0[r][i] for r in range(k)) != Fraction(1):
            raise CertificateError("B0 is not doubly stochastic")
        if any(B0[i][j] <= 0 for j in range(k)):
            raise CertificateError("B0 must be strictly positive (increase M_param)")
    if mat_det_exact(B0) == 0:
        raise CertificateError("B0 must be invertible")

    # Verify exact projective bijection y = h*x / (1 + c^T x), inverse x = y / (h - c^T y),
    # and exact slack identity 1 - B0*y == (1 - A*x) / (1 + c^T x) on test_x!
    raw_x = payload.get("test_x")
    if not isinstance(raw_x, list) or len(raw_x) != k:
        raise CertificateError(f"test_x must have length {k}")
    x_vec = [_parse_frac_pair(v, f"test_x[{j}]") for j, v in enumerate(raw_x)]
    if any(xv < 0 for xv in x_vec) or all(xv == 0 for xv in x_vec):
        raise CertificateError("test_x must be a nonzero nonnegative vector")
    Ax = [sum(A[i][j] * x_vec[j] for j in range(k)) for i in range(k)]
    if any(v > Fraction(1) for v in Ax):
        raise CertificateError("test_x is infeasible for A*x <= 1")
    cx = sum(c_vec[j] * x_vec[j] for j in range(k))
    denom_x = Fraction(1) + cx
    y_vec = [frac_div(h * xv, denom_x) for xv in x_vec]
    cy = sum(c_vec[j] * y_vec[j] for j in range(k))
    if cy >= h:
        raise CertificateError("Projective image violated c^T y < h")
    x_back = [frac_div(yv, h - cy) for yv in y_vec]
    if x_back != x_vec:
        raise CertificateError("Stage-1 projective inverse failed to recover test_x")
    B0y = [sum(B0[i][j] * y_vec[j] for j in range(k)) for i in range(k)]
    for i in range(k):
        if Fraction(1) - B0y[i] != frac_div(Fraction(1) - Ax[i], denom_x):
            raise CertificateError("Stage-1 projective slack identity 1 - B0*y = (1 - A*x)/(1 + c^T x) failed")

    # Stage 2: Symmetric lifting + projective smoothing B = (H + t*1*1^T)/(1 + d*t)
    t_sm = _req_frac(payload, "t_smooth")
    eps = _req_frac(payload, "epsilon_bound")
    if t_sm <= 0 or eps <= 0:
        raise CertificateError("t_smooth and epsilon_bound must be positive")
    d = 2 * k + (1 if odd_pad else 0)
    if _req_int(payload, "claimed_d") != d:
        raise CertificateError(f"claimed_d != exact {d}")

    B = stage2_symmetric_smooth_packing(B0, odd_pad, t_sm)
    target = Fraction(1, d)
    max_dev = Fraction(0)
    for i in range(d):
        if sum(B[i]) != Fraction(1):
            raise CertificateError(f"Stage-2 B row {i} sum != 1")
        for j in range(d):
            if B[i][j] != B[j][i]:
                raise CertificateError("Stage-2 B is not symmetric")
            if B[i][j] <= 0:
                raise CertificateError("Stage-2 B is not strictly positive")
            dev = abs(B[i][j] - target)
            if dev > max_dev:
                max_dev = dev
    if mat_det_exact(B) == 0:
        raise CertificateError("Stage-2 symmetric doubly stochastic B must be invertible")
    if max_dev >= eps:
        raise CertificateError(f"Stage-2 max deviation {max_dev} >= epsilon_bound {eps}")
    if _req_frac(payload, "claimed_max_dev") != max_dev:
        raise CertificateError(f"claimed_max_dev disagrees with exact {max_dev}")

    return {
        "mode": "symmetric_doubly_stochastic_packing_cert",
        "k_dim": k,
        "final_d": d,
        "max_dev_from_1_over_d": f"{max_dev.numerator}/{max_dev.denominator}",
    }


# ---------------------------------------------------------------------------
# Public Entry Points
# ---------------------------------------------------------------------------

def verify_certificate_dict(data: Mapping[str, Any]) -> Dict[str, Any]:
    if not isinstance(data, Mapping):
        raise CertificateError("Certificate must be a mapping")
    allowed_top = {"protocol", "certificate_id", "mode", "payload", "description", "note"}
    extra_top = set(data.keys()) - allowed_top
    if extra_top:
        raise CertificateError(f"Unexpected top-level keys: {sorted(extra_top)}")
    for opt_k in ("description", "note"):
        if opt_k in data and not isinstance(data[opt_k], str):
            raise CertificateError(f"Optional field {opt_k!r} must be a string")

    proto = _req_str(data, "protocol")
    if proto != PROTOCOL_ID:
        raise CertificateError(f"Unsupported protocol {proto!r}; expected {PROTOCOL_ID!r}")
    cert_id = _req_str(data, "certificate_id")
    mode = _req_str(data, "mode")
    if mode not in VALID_MODES:
        raise CertificateError(f"Unknown certificate mode {mode!r}")

    payload = data.get("payload")
    if not isinstance(payload, Mapping):
        raise CertificateError("Field 'payload' must be a JSON object")

    try:
        if mode == "convex_corner_bend_and_crossing_cert":
            summary = _verify_mode1(payload)
        elif mode == "klee_walkup_wedge_and_complement_cert":
            summary = _verify_mode2(payload)
        elif mode == "schlegel_projective_and_prismatoid_cert":
            summary = _verify_mode3(payload)
        else:
            summary = _verify_mode4(payload)
    except CertificateError:
        raise
    except ValueError as exc:
        raise CertificateError(f"Exact verification value error: {exc}") from exc

    summary["certificate_id"] = cert_id
    summary["protocol"] = PROTOCOL_ID
    return summary


def verify_certificate_str(raw: str) -> Dict[str, Any]:
    return verify_certificate_dict(load_strict_json_str(raw))


def verify_certificate_file(path: Path | str) -> Dict[str, Any]:
    return verify_certificate_str(Path(path).read_text(encoding="utf-8"))
