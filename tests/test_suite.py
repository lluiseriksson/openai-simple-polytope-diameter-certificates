"""Comprehensive 82-test suite for HIRSCH-DIAM-CERT-1.0.

Includes:
  - 1 aggregate canonical 16-certificate replay,
  - 3 raw-JSON zero-float / NaN / duplicate-key parser tests,
  - 2 exact arithmetic soundness regressions (Cor. 14.3 Stage-2 symmetry/invertibility & 4D hexagon wedge counts),
  - 76 parametrized live adversarial forgery tests (12 schema/parser + 15 Mode 1 + 16 Mode 2 +
    18 Mode 3 + 15 Mode 4), each asserting an exact regex match on the intended CertificateError branch.
"""

from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
import pytest

from src.hirsch_certifier import (
    PROTOCOL_ID,
    verify_certificate_dict,
    verify_certificate_file,
    verify_certificate_str,
)
from src.hirsch_exact import (
    CertificateError,
    mat_det_exact,
    stage1_doubly_stochastic,
    stage2_symmetric_smooth_packing,
)

ROOT = Path(__file__).resolve().parent.parent
CERT_DIR = ROOT / "certificates"


def _load(name: str) -> dict:
    return json.loads((CERT_DIR / f"{name}.cert.json").read_text(encoding="utf-8"))


def test_all_16_canonical_certificates_verify() -> None:
    files = sorted(CERT_DIR.glob("*.cert.json"))
    assert len(files) == 16
    for f in files:
        res = verify_certificate_file(f)
        assert res["protocol"] == PROTOCOL_ID


def test_raw_json_rejects_float_literal() -> None:
    raw = '{"protocol": "HIRSCH-DIAM-CERT-1.0", "certificate_id": "x", "mode": "convex_corner_bend_and_crossing_cert", "payload": {"N": 16.5}}'
    with pytest.raises(CertificateError, match="Floating-point literal rejected"):
        verify_certificate_str(raw)


def test_raw_json_rejects_nan_constant() -> None:
    raw = '{"protocol": "HIRSCH-DIAM-CERT-1.0", "certificate_id": "x", "mode": "convex_corner_bend_and_crossing_cert", "payload": {"N": NaN}}'
    with pytest.raises(CertificateError, match="Non-standard JSON constant rejected"):
        verify_certificate_str(raw)


def test_raw_json_rejects_duplicate_keys() -> None:
    raw = '{"protocol": "HIRSCH-DIAM-CERT-1.0", "protocol": "HIRSCH-DIAM-CERT-1.0"}'
    with pytest.raises(CertificateError, match="Duplicate JSON key rejected"):
        verify_certificate_str(raw)


def test_regression_cor143_two_stage_invertible_and_symmetric() -> None:
    A = [[Fraction(3, 4), Fraction(1, 4)], [Fraction(1, 5), Fraction(4, 5)]]
    _, _, B0 = stage1_doubly_stochastic(A, Fraction(2))
    B = stage2_symmetric_smooth_packing(B0, True, Fraction(10))
    assert len(B) == 5
    assert mat_det_exact(B) != 0
    for i in range(5):
        assert sum(B[i]) == Fraction(1)
        for j in range(5):
            assert B[i][j] == B[j][i] and B[i][j] > 0


def test_regression_wedge_complementary_count() -> None:
    res = verify_certificate_file(CERT_DIR / "m2_wedge_hexagon_2d_to_4d.cert.json")
    assert res["final_dim"] == 4 and res["final_facets"] == 8


def _build_76_forgeries() -> list[tuple[str, dict, str]]:
    m1 = _load("m1_corners_3_n16_k2")
    m2 = _load("m2_wedge_pentagon_2d_to_3d")
    m3 = _load("m3_schlegel_antiprism_3d_fig7")
    m3_4d = _load("m3_schlegel_antiprism_4d")
    m4 = _load("m4_packing_k2_even_d4")

    cases: list[tuple[str, dict, str]] = []

    def add(name: str, base: dict, mutator, regex: str) -> None:
        c = copy.deepcopy(base)
        mutator(c)
        cases.append((name, c, regex))

    # --- Schema / Top-Level (12) ---
    add("S01_bad_protocol", m1, lambda c: c.__setitem__("protocol", "FAKE-1.0"), "Unsupported protocol")
    add("S02_whitespace_cert_id", m1, lambda c: c.__setitem__("certificate_id", "  "), "non-empty string")
    add("S03_unknown_mode", m1, lambda c: c.__setitem__("mode", "bogus_mode"), "Unknown certificate mode")
    add("S04_extra_top_key", m1, lambda c: c.__setitem__("rogue", 1), "Unexpected top-level keys")
    add("S05_bad_description", m1, lambda c: c.__setitem__("description", 99), "Optional field 'description' must be a string")
    add("S06_bad_note", m1, lambda c: c.__setitem__("note", True), "Optional field 'note' must be a string")
    add("S07_payload_not_obj", m1, lambda c: c.__setitem__("payload", "bad"), "must be a JSON object")
    add("S08_missing_protocol", m1, lambda c: c.pop("protocol"), "Missing required string field 'protocol'")
    add("S09_bool_as_int", m1, lambda c: c["payload"].__setitem__("N_param", True), "must be an exact integer")
    add("S10_float_in_dict", m1, lambda c: c["payload"].__setitem__("K_param", 2.0), "must be an exact integer")
    add("S11_missing_key", m1, lambda c: c["payload"].pop("corner_xs"), "Missing keys in Mode 1 payload")
    add("S12_extra_key", m1, lambda c: c["payload"].__setitem__("extra", 1), "Unexpected keys in Mode 1 payload")

    # --- Mode 1 Forgeries (15) ---
    add("M1_01_single_corner", m1, lambda c: c["payload"].__setitem__("corner_xs", [[1, 1]]), "2\\.\\.12 strictly increasing")
    add("M1_02_corners_unsorted", m1, lambda c: c["payload"]["corner_xs"].reverse(), "strictly increasing")
    add("M1_03_unreduced_frac", m1, lambda c: c["payload"]["corner_xs"].__setitem__(0, [2, 2]), "lowest terms \\(gcd == 1\\)")
    add("M1_04_path_not_spanning", m1, lambda c: c["payload"]["path_points"][0].__setitem__(0, [1, 1]), "x-coordinates must be strictly increasing")
    add("M1_05_corner_nonzero_z", m1, lambda c: c["payload"]["path_points"][1].__setitem__(1, [1, 10]), "must attain z=0")
    add("M1_06_corner_flat_slope", m1, lambda c: c["payload"]["path_points"][0].__setitem__(1, [1, 10]), "violates V-corner lower envelope")
    add("M1_07_wrong_breaks", m1, lambda c: c["payload"].__setitem__("claimed_breaks", 3), "claimed_breaks=3 disagrees")
    add("M1_08_N_too_small", m1, lambda c: c["payload"].__setitem__("N_param", 4), "parameters \\(N, K, R\\) out")
    add("M1_09_K_too_large", m1, lambda c: c["payload"].__setitem__("K_param", 9), "parameters \\(N, K, R\\) out")
    add("M1_10_gap_gt_one", m1, lambda c: c["payload"].__setitem__("gap_g", [3, 2]), "gap_g must lie in")
    add("M1_11_wrong_full_lb", m1, lambda c: c["payload"].__setitem__("claimed_full_dist_lb", 99), "claimed_full_dist_lb disagrees")
    add("M1_12_wrong_gap_lb", m1, lambda c: c["payload"].__setitem__("claimed_gap_dist_lb", [1, 1]), "claimed_gap_dist_lb disagrees")
    add("M1_13_bad_prod_d", m1, lambda c: c["payload"].__setitem__("prod_d", 4), "product padding dimensions/facets invalid")
    add("M1_14_bad_prod_E", m1, lambda c: c["payload"].__setitem__("prod_E", 1), "product padding dimensions/facets invalid")
    add("M1_15_neg_z_waypoint", m1, lambda c: c["payload"]["path_points"].insert(2, [[3, 2], [-1000, 1]]), "violates V-corner lower envelope")

    # --- Mode 2 Forgeries (16) ---
    add("M2_01_bad_base_dim", m2, lambda c: c["payload"].__setitem__("base_dim", 4), "base_dim must be 2 or 3")
    add("M2_02_rhs_len_mismatch", m2, lambda c: c["payload"]["base_rhs"].pop(), "lists of equal length")
    add("M2_03_too_few_facets", m2, lambda c: (c["payload"]["base_normals"].pop(), c["payload"]["base_rhs"].pop()), "uncovered facets 1\\.\\.3")
    add("M2_04_degenerate_polytope", m2, lambda c: c["payload"]["base_rhs"].__setitem__(4, [2, 1]), "Degenerate vertex")
    add("M2_05_same_tips", m2, lambda c: c["payload"].__setitem__("tip_w_idx", c["payload"]["tip_v_idx"]), "distinct valid base vertex indices")
    add("M2_06_tip_oob", m2, lambda c: c["payload"].__setitem__("tip_w_idx", 99), "distinct valid base vertex indices")
    add("M2_07_adjacent_tips", m2, lambda c: c["payload"].__setitem__("tip_w_idx", 1), "share a facet")
    add("M2_08_wedge_seq_len", m2, lambda c: c["payload"].__setitem__("wedge_facet_sequence", [3, 0]), "uncovered facets 1")
    add("M2_09_wedge_covered_facet", m2, lambda c: c["payload"].__setitem__("wedge_facet_sequence", [0]), "not strictly uncovered by both tips")
    add("M2_10_wedge_idx_oob", m2, lambda c: c["payload"].__setitem__("wedge_facet_sequence", [9]), "out of range")
    add("M2_11_wrong_final_dim", m2, lambda c: c["payload"].__setitem__("claimed_final_dim", 4), "claimed_final_dim mismatch")
    add("M2_12_wrong_final_facets", m2, lambda c: c["payload"].__setitem__("claimed_final_facets", 7), "claimed_final_facets mismatch")
    add("M2_13_wrong_final_verts", m2, lambda c: c["payload"].__setitem__("claimed_final_vertex_count", 99), "claimed_final_vertex_count")
    add("M2_14_unreduced_normal", m2, lambda c: c["payload"]["base_normals"][0].__setitem__(0, [-2, 2]), "lowest terms")
    add(
        "M2_15_zero_normal_pseudofacet",
        m2,
        lambda c: c["payload"]["base_normals"].__setitem__(3, [[0, 1], [0, 1]]),
        "cannot be the zero vector",
    )
    add(
        "M2_16_unbounded_polyhedron",
        m2,
        lambda c: (
            c["payload"].__setitem__(
                "base_normals",
                [[[-1, 1], [0, 1]], [[0, 1], [-1, 1]], [[-2, 1], [-1, 1]], [[-1, 1], [-2, 1]], [[-1, 1], [-1, 1]]],
            ),
            c["payload"].__setitem__("base_rhs", [[0, 1], [0, 1], [-2, 1], [-2, 1], [10, 1]]),
        ),
        "unbounded or non-polytopal",
    )

    # --- Mode 3 Forgeries (18) ---
    add("M3_01_bad_dm1", m3, lambda c: c["payload"].__setitem__("base_dim_d_minus_1", 1), "base_dim_d_minus_1 must lie in")
    add("M3_02_T0_wrong_count", m3, lambda c: c["payload"]["T0_vertices"].pop(), "exactly d=3 rational simplex vertices")
    add("M3_03_T0_dup_vert", m3, lambda c: c["payload"]["T0_vertices"].__setitem__(2, list(c["payload"]["T0_vertices"][0])), "vertices must be distinct")
    add("M3_04_T0_collinear", m3, lambda c: c["payload"]["T0_vertices"].__setitem__(2, [[-1, 1], [2, 1]]), "nondegenerate 2-simplex")
    add("M3_05_wrong_minkowski_cnt", m3, lambda c: c["payload"].__setitem__("claimed_minkowski_vertex_count", 3), "claimed_minkowski_vertex_count")
    add("M3_06_empty_facets", m3, lambda c: c["payload"].__setitem__("nonbase_facets", []), "3\\.\\.32 non-base inequalities")
    add("M3_07_facet_bad_key", m3, lambda c: c["payload"]["nonbase_facets"][0].__setitem__("d", [0, 1]), "Unexpected keys in nonbase_facets\\[0\\]")
    add("M3_08_viewpoint_not_strict", m3, lambda c: c["payload"]["nonbase_facets"][0].__setitem__("b", [-1, 1]), "Schlegel viewpoint strictness")
    add("M3_09_facet_infeasible_T0", m3, lambda c: c["payload"]["nonbase_facets"][0].__setitem__("c", [1, 2]), "violated by T0 vertex")
    add("M3_10_facet_not_supporting", m3, lambda c: c["payload"]["nonbase_facets"][0].__setitem__("c", [5, 1]), "not a supporting \\(d-1\\)-facet")
    add("M3_11_Ds_inflated_cheat", m3, lambda c: (c["payload"].__setitem__("prismatoid_width_Ds", 99), c["payload"].__setitem__("claimed_section_dist", 97), c["payload"].__setitem__("claimed_schlegel_depth", 98)), "disagrees with exact dual BFS width 3")
    add("M3_12_wrong_section_dist", m3, lambda c: c["payload"].__setitem__("claimed_section_dist", 2), "claimed_section_dist")
    add("M3_13_wrong_schlegel_depth", m3, lambda c: c["payload"].__setitem__("claimed_schlegel_depth", 3), "claimed_schlegel_depth")
    add("M3_14_facet_infeasible_T1", m3, lambda c: c["payload"]["nonbase_facets"][0].__setitem__("b", [1, 1]), "violated by T1 vertex")
    add("M3_15_base_facet_T1_cheat", m3, lambda c: c["payload"]["nonbase_facets"].__setitem__(0, {"a": [[0, 1], [0, 1]], "b": [1, 1], "c": [1, 1]}), "zero horizontal normal")
    add("M3_16_duplicate_facet", m3, lambda c: c["payload"]["nonbase_facets"].__setitem__(1, copy.deepcopy(c["payload"]["nonbase_facets"][0])), "Duplicate non-base facet rejected")
    add("M3_17_incomplete_facet_list", m3_4d, lambda c: c["payload"]["nonbase_facets"].pop(), "Incomplete facet list")
    add("M3_18_origin_outside_T0", m3, lambda c: c["payload"].__setitem__("T0_vertices", [[[1, 1], [1, 1]], [[3, 1], [1, 1]], [[1, 1], [3, 1]]]), "Origin 0 must lie strictly in the interior")

    # --- Mode 4 Forgeries (15) ---
    add("M4_01_k_too_small", m4, lambda c: c["payload"].__setitem__("k_dim", 1), "k_dim must lie in")
    add("M4_02_odd_pad_not_bool", m4, lambda c: c["payload"].__setitem__("odd_pad", 0), "odd_pad must be a boolean")
    add("M4_03_A_row_sum_not_one", m4, lambda c: c["payload"]["A_matrix"][0].__setitem__(0, [1, 2]), "violates A\\*1 = 1")
    add("M4_04_A_singular", m4, lambda c: c["payload"]["A_matrix"].__setitem__(1, list(c["payload"]["A_matrix"][0])), "must be invertible")
    add("M4_05_M_nonpositive", m4, lambda c: c["payload"].__setitem__("M_param", [0, 1]), "M_param must be positive")
    add("M4_06_M_too_small_neg_c", m4, lambda c: c["payload"].__setitem__("M_param", [1, 4]), "Stage-1 vector c must be strictly positive")
    add("M4_07_test_x_all_zero", m4, lambda c: c["payload"].__setitem__("test_x", [[0, 1], [0, 1]]), "nonzero nonnegative vector")
    add("M4_08_test_x_infeasible", m4, lambda c: c["payload"].__setitem__("test_x", [[2, 1], [2, 1]]), "infeasible for A\\*x <= 1")
    add("M4_09_test_x_neg", m4, lambda c: c["payload"]["test_x"].__setitem__(0, [-1, 2]), "nonzero nonnegative vector")
    add("M4_10_t_smooth_zero", m4, lambda c: c["payload"].__setitem__("t_smooth", [0, 1]), "t_smooth and epsilon_bound must be positive")
    add("M4_11_wrong_claimed_d", m4, lambda c: c["payload"].__setitem__("claimed_d", 5), "claimed_d != exact 4")
    add("M4_12_eps_too_tight", m4, lambda c: c["payload"].__setitem__("epsilon_bound", [1, 1000]), "Stage-2 max deviation")
    add("M4_13_wrong_max_dev", m4, lambda c: c["payload"].__setitem__("claimed_max_dev", [1, 100]), "claimed_max_dev disagrees")
    add("M4_14_unreduced_t_smooth", m4, lambda c: c["payload"].__setitem__("t_smooth", [10, 2]), "lowest terms")
    add("M4_15_A_negative_nonpacking", m4, lambda c: c["payload"].__setitem__("A_matrix", [[[2, 1], [-1, 1]], [[-1, 1], [2, 1]]]), "violates packing LP A >= 0")

    assert len(cases) == 76
    return cases


FORGERIES = _build_76_forgeries()


@pytest.mark.parametrize("name,cert,expected_regex", FORGERIES, ids=[c[0] for c in FORGERIES])
def test_adversarial_forgery_rejected(name: str, cert: dict, expected_regex: str) -> None:
    with pytest.raises(CertificateError, match=expected_regex):
        verify_certificate_dict(cert)
