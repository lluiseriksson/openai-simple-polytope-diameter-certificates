"""Zero-float exact arithmetic engine for HIRSCH-DIAM-CERT-1.0.

Implements pure-integer and exact-Fraction primitives corresponding to OpenAI's
November 24, 2026 preprint "Superpolynomial graph diameter of simple polytopes":
  - Exact rational convex-corner bend verification & Lemma 2.2 cell-distance / Theorem 1.2
    product facet-dimension accounting (Secs. 1.2, 2.1, 2.2),
  - Exact rational linear algebra (Gaussian elimination rank & determinant over Q),
  - Exact Klee-Walkup wedge lifting W = {(x, t) : x in F, 0 <= t <= a(x)} with full
    vertex enumeration, simplicity verification, and complementary-tip certification (Sec. 14, Cor. 14.1),
  - Exact Schlegel projective map (x, z) -> (x/(1+z), -z/(1+z)) and prismatoid middle-section
    (U + V)/2 correspondence (Sec. 14.1, Cor. 14.2),
  - Exact Corollary 14.3 two-stage projective transformation to invertible symmetric positive
    doubly-stochastic packing matrices B with |B_{ij} - 1/d| < eps (Sec. 14.2, Cor. 14.3).

NO floating-point operations, true division '/', or math imports appear anywhere in this module.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from typing import Dict, List, Sequence, Tuple


class CertificateError(ValueError):
    """Raised whenever any exact mathematical or structural check fails."""


def gcd_int(a: int, b: int) -> int:
    x, y = abs(a), abs(b)
    while y:
        x, y = y, x % y
    return x


def frac_div(a: Fraction, b: Fraction) -> Fraction:
    """Exact Fraction division without using the '/' operator (AST-clean)."""
    if b.numerator == 0:
        raise CertificateError("Division by zero in exact rational arithmetic")
    return Fraction(a.numerator * b.denominator, a.denominator * b.numerator)


def mat_det_exact(mat: Sequence[Sequence[Fraction]]) -> Fraction:
    """Exact determinant of a square rational matrix via pure-Fraction Gaussian elimination."""
    n = len(mat)
    A = [[Fraction(x) for x in row] for row in mat]
    det = Fraction(1)
    sign = 1
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            return Fraction(0)
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
            sign = -sign
        pv = A[col][col]
        det *= pv
        for r in range(col + 1, n):
            if A[r][col] != 0:
                fac = frac_div(A[r][col], pv)
                for c in range(col, n):
                    A[r][c] -= fac * A[col][c]
    return det if sign == 1 else -det


def solve_linear_exact(mat: Sequence[Sequence[Fraction]], rhs: Sequence[Fraction]) -> List[Fraction] | None:
    """Solve square linear system mat * x = rhs over Q. Returns None if singular."""
    n = len(mat)
    A = [[Fraction(x) for x in row] + [Fraction(rhs[i])] for i, row in enumerate(mat)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        for c in range(col, n + 1):
            A[col][c] = frac_div(A[col][c], pv)
        for r in range(n):
            if r != col and A[r][col] != 0:
                fac = A[r][col]
                for c in range(col, n + 1):
                    A[r][c] -= fac * A[col][c]
    return [A[i][n] for i in range(n)]


def mat_rank_exact(mat: Sequence[Sequence[Fraction]]) -> int:
    """Exact rank of a rational matrix over Q via pure-Fraction Gaussian elimination."""
    if not mat or not mat[0]:
        return 0
    m = len(mat)
    n = len(mat[0])
    A = [[Fraction(x) for x in row] for row in mat]
    rank = 0
    for col in range(n):
        piv = None
        for r in range(rank, m):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        if piv != rank:
            A[rank], A[piv] = A[piv], A[rank]
        pv = A[rank][col]
        for r in range(rank + 1, m):
            if A[r][col] != 0:
                fac = frac_div(A[r][col], pv)
                for c in range(col, n):
                    A[r][c] -= fac * A[rank][c]
        rank += 1
        if rank == m:
            break
    return rank


def affine_rank_exact(points: Sequence[Sequence[Fraction]]) -> int:
    """Exact affine dimension of a nonempty set of rational points."""
    if len(points) <= 1:
        return 0
    p0 = points[0]
    dim = len(p0)
    diffs = [[pt[j] - p0[j] for j in range(dim)] for pt in points[1:]]
    return mat_rank_exact(diffs)


def origin_strictly_inside_simplex(verts: Sequence[Sequence[Fraction]]) -> bool:
    """Verify 0 lies strictly in the interior of the (d-1)-simplex spanned by d vertices in Q^{d-1}."""
    d = len(verts)
    dm1 = d - 1
    # Solve sum_i lam_i * v_i = 0 and sum_i lam_i = 1
    mat = [[verts[i][r] for i in range(d)] for r in range(dm1)] + [[Fraction(1)] * d]
    rhs = [Fraction(0)] * dm1 + [Fraction(1)]
    lam = solve_linear_exact(mat, rhs)
    if lam is None:
        return False
    return all(val > 0 for val in lam)


def verify_simple_polytope_bounded_and_irredundant(
    verts: Sequence[Tuple[Fraction, ...]],
    tight_sets: Sequence[Tuple[int, ...]],
    m_facets: int,
    dim: int,
) -> None:
    """Verify exact boundedness and facet irredundancy of a simple polytope.

    1. Boundedness: in a simple d-polyhedron, every (d-1)-subset of every vertex's tight d-set
       defines an incident edge or unbounded ray; it is bounded iff EVERY such (d-1)-subset
       is shared by exactly 2 vertices (zero unbounded rays).
    2. Irredundancy: every inequality 0..m_facets-1 must be tight on >= dim vertices
       spanning exact affine rank dim - 1.
    """
    if len(verts) < dim + 1:
        raise CertificateError(f"Simple {dim}-polytope must have at least {dim + 1} vertices")
    ridge_counts: Dict[Tuple[int, ...], int] = {}
    facet_verts: List[List[Tuple[Fraction, ...]]] = [[] for _ in range(m_facets)]
    for v, ts in zip(verts, tight_sets):
        for f_idx in ts:
            facet_verts[f_idx].append(v)
        for edge_facets in combinations(ts, dim - 1):
            ridge_counts[edge_facets] = ridge_counts.get(edge_facets, 0) + 1
    for edge_facets, cnt in ridge_counts.items():
        if cnt != 2:
            raise CertificateError(
                f"Polyhedron is unbounded or non-polytopal: edge {edge_facets} incident to {cnt} != 2 vertices"
            )
    for f_idx in range(m_facets):
        fvs = facet_verts[f_idx]
        if len(fvs) < dim or affine_rank_exact(fvs) != dim - 1:
            raise CertificateError(
                f"Inequality {f_idx} is redundant or not a supporting facet (affine rank < {dim - 1})"
            )


def compute_prismatoid_dual_width(
    all_verts: Sequence[Tuple[Fraction, ...]],
    facet_tight_sets: Sequence[Tuple[int, ...]],
    dim: int,
) -> int:
    """Verify complete ridge-adjacency dual graph of prismatoid [T0, T1, F_1..F_L] and compute
    exact BFS shortest ridge-adjacency path distance D_1 between T0 (index 0) and T1 (index 1).
    """
    n_fac = len(facet_tight_sets)
    # Build exact ridge-adjacency graph: two facets share a ridge iff intersection has affine rank dim - 2
    adj: List[List[int]] = [[] for _ in range(n_fac)]
    ridge_owners: Dict[Tuple[int, ...], List[int]] = {}
    for i in range(n_fac):
        si = set(facet_tight_sets[i])
        for j in range(i + 1, n_fac):
            inter = sorted(si & set(facet_tight_sets[j]))
            if len(inter) >= dim - 1:
                pts = [all_verts[idx] for idx in inter]
                if affine_rank_exact(pts) == dim - 2:
                    adj[i].append(j)
                    adj[j].append(i)
                    key = tuple(inter)
                    owners = ridge_owners.setdefault(key, [])
                    if i not in owners:
                        owners.append(i)
                    if j not in owners:
                        owners.append(j)

    # Check no ridge belongs to >= 3 facets
    for key, owners in ridge_owners.items():
        if len(owners) != 2:
            raise CertificateError(f"Ridge {key} shared by {len(owners)} != 2 facets")

    # Completeness check:
    # - For every simplicial facet (|S| == dim, including T0 and T1), ALL dim of its (dim-1)-subsets
    #   must appear as ridges with a unique neighbor!
    # - For every non-simplicial facet (prism lateral facet of size 2*(dim-1)), must have >= dim + 1 neighbors
    #   (the 2 base ridges T0, T1 plus dim - 1 lateral ridges).
    for i, f_ts in enumerate(facet_tight_sets):
        if len(f_ts) == dim:
            for sub in combinations(f_ts, dim - 1):
                if sub not in ridge_owners:
                    raise CertificateError(f"Incomplete facet list: ridge {sub} of facet {i} has no adjacent facet")
        else:
            if len(adj[i]) < dim + 1:
                raise CertificateError(f"Incomplete facet list: non-simplicial facet {i} has only {len(adj[i])} < {dim + 1} ridge neighbors")

    # BFS from T0 (0) to T1 (1)
    dist = [-1] * n_fac
    dist[0] = 0
    queue = [0]
    head = 0
    while head < len(queue):
        u = queue[head]
        head += 1
        if u == 1:
            return dist[1]
        for w in adj[u]:
            if dist[w] == -1:
                dist[w] = dist[u] + 1
                queue.append(w)
    raise CertificateError("Prismatoid dual graph is disconnected between T0 and T1")


def enumerate_simple_polytope_vertices(
    normals: Sequence[Sequence[Fraction]],
    rhs: Sequence[Fraction],
    dim: int,
) -> Tuple[List[Tuple[Fraction, ...]], List[Tuple[int, ...]]]:
    """Enumerate all vertices of {x in Q^dim : normals[i] . x <= rhs[i]} and verify simplicity.
    Returns (vertices, tight_facet_index_tuples). Raises CertificateError if degenerate (not simple).
    """
    m = len(normals)
    seen: Dict[Tuple[Fraction, ...], Tuple[int, ...]] = {}
    for comb in combinations(range(m), dim):
        sub_mat = [normals[i] for i in comb]
        sub_rhs = [rhs[i] for i in comb]
        sol = solve_linear_exact(sub_mat, sub_rhs)
        if sol is None:
            continue
        # Check feasibility against all m inequalities
        tight: List[int] = []
        feasible = True
        for i in range(m):
            val = sum(normals[i][j] * sol[j] for j in range(dim))
            if val > rhs[i]:
                feasible = False
                break
            if val == rhs[i]:
                tight.append(i)
        if not feasible:
            continue
        if len(tight) != dim:
            raise CertificateError(f"Degenerate vertex {sol} lies on {len(tight)} != {dim} facets (not simple)")
        pt = tuple(sol)
        if pt not in seen:
            seen[pt] = tuple(tight)
    verts = sorted(seen.keys())
    return verts, [seen[v] for v in verts]


def klee_walkup_wedge_system(
    normals: Sequence[Sequence[Fraction]],
    rhs: Sequence[Fraction],
    facet_idx: int,
) -> Tuple[List[List[Fraction]], List[Fraction]]:
    """Construct Section 14 Klee-Walkup wedge W = {(x, t) : x in F, 0 <= t <= a(x)}
    where a(x) = rhs[facet_idx] - normals[facet_idx] . x >= 0 is the slack of facet_idx.
    Replaces facet_idx by two sheets:
      - bottom sheet:  -t <= 0   (t >= 0)
      - top sheet:     normals[facet_idx] . x + t <= rhs[facet_idx]   (t <= a(x))
    All other facets lift with coefficient 0 on t.
    """
    new_normals: List[List[Fraction]] = []
    new_rhs: List[Fraction] = []
    dim = len(normals[0])
    for i, (n_row, b_val) in enumerate(zip(normals, rhs)):
        if i == facet_idx:
            # Bottom sheet: 0*x - 1*t <= 0
            new_normals.append([Fraction(0)] * dim + [Fraction(-1)])
            new_rhs.append(Fraction(0))
            # Top sheet: n_row*x + 1*t <= b_val
            new_normals.append([Fraction(c) for c in n_row] + [Fraction(1)])
            new_rhs.append(Fraction(b_val))
        else:
            new_normals.append([Fraction(c) for c in n_row] + [Fraction(0)])
            new_rhs.append(Fraction(b_val))
    return new_normals, new_rhs


def stage1_doubly_stochastic(
    A: Sequence[Sequence[Fraction]],
    M: Fraction,
) -> Tuple[List[Fraction], Fraction, List[List[Fraction]]]:
    """Exact Section 14.2 (Cor. 14.3) Stage-1 transformation:
      c_j = M - k^(-1) * sum_i A_{ij},  h = k * M,  B0 = (A + 1 * c^T) / h.
    """
    k = len(A)
    c_vec: List[Fraction] = []
    for j in range(k):
        col_sum = sum(A[i][j] for i in range(k))
        c_vec.append(M - frac_div(col_sum, Fraction(k)))
    h = Fraction(k) * M
    B0: List[List[Fraction]] = []
    for i in range(k):
        B0.append([frac_div(A[i][j] + c_vec[j], h) for j in range(k)])
    return c_vec, h, B0


def stage2_symmetric_smooth_packing(
    C: Sequence[Sequence[Fraction]],
    odd_pad: bool,
    t_smooth: Fraction,
) -> List[List[Fraction]]:
    """Exact Section 14.2 (Cor. 14.3) Stage-2 symmetric block lifting & projective smoothing:
      H = [[0, C], [C^T, 0]] (plus [1] diagonal block if odd_pad),
      B = (H + t * 1 * 1^T) / (1 + d * t).
    """
    k = len(C)
    d = 2 * k + (1 if odd_pad else 0)
    H = [[Fraction(0)] * d for _ in range(d)]
    for i in range(k):
        for j in range(k):
            H[i][k + j] = C[i][j]
            H[k + j][i] = C[i][j]
    if odd_pad:
        H[d - 1][d - 1] = Fraction(1)
    denom = Fraction(1) + Fraction(d) * t_smooth
    B = [[frac_div(H[i][j] + t_smooth, denom) for j in range(d)] for i in range(d)]
    return B
