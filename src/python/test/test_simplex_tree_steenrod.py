""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Anibal M. Medina-Mardones

    Copyright (C) 2026 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import math
from itertools import combinations

import pytest

from gudhi import SimplexTree


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _all_faces(tops):
    """Return all faces (as sorted tuples) of a list of top simplices."""
    faces = set()
    for s in tops:
        for d in range(len(s) + 1):
            for f in combinations(s, d + 1):
                faces.add(tuple(sorted(f)))
    return faces


def _simplex_tree_with_integer_filtration(faces, offset=0):
    """Insert *faces* into a SimplexTree using contiguous integer filtration
    values starting at *offset*.

    Simplices are ordered by (dimension, lex) so that every face of a simplex
    receives a strictly smaller filtration value, giving a valid filtration.
    The filtration value of the i-th simplex (0-indexed) is ``offset + i``.

    Returns ``(simplex_tree, n_simplices)``.
    """
    ordered = sorted(faces, key=lambda f: (len(f), f))
    st = SimplexTree()
    for i, f in enumerate(ordered):
        st.insert(list(f), filtration=float(offset + i))
    return st, len(ordered)


# ---------------------------------------------------------------------------
# RP² — Sq¹, essential bars
# ---------------------------------------------------------------------------

# Minimal triangulation of RP²: 6 vertices, 15 edges, 10 triangles.
# H*(RP²; F_2) = F_2 in degrees 0, 1, 2.  Sq^1: H^1 -> H^2 is an isomorphism.
RP2_TOP = [
    [1, 2, 4], [2, 3, 4], [1, 3, 5], [2, 3, 5], [1, 4, 5],
    [1, 2, 6], [1, 3, 6], [3, 4, 6], [2, 5, 6], [4, 5, 6],
]


def _rp2_simplex_tree():
    st, _ = _simplex_tree_with_integer_filtration(_all_faces(RP2_TOP))
    return st


def test_rp2_sq1_three_essential_ordinary_bars():
    """RP² has three essential ordinary bars (H^0, H^1, H^2 over F_2)."""
    ordinary, _ = _rp2_simplex_tree().compute_steenrod_barcodes(k=1)
    n_ess = sum(1 for d in ordinary for _, death in d if math.isinf(death))
    assert n_ess == 3


def test_rp2_sq1_one_essential_steenrod_bar_in_dim2():
    """Sq^1 produces exactly one essential bar, in cohomological degree 2.

    This reflects the non-orientability of RP^2: Sq^1: H^1 -> H^2 is an
    isomorphism over F_2.
    """
    _, steenrod = _rp2_simplex_tree().compute_steenrod_barcodes(k=1)
    sq1_ess = [b for b in steenrod[2] if math.isinf(b[1])]
    assert len(sq1_ess) == 1


def test_does_not_require_compute_persistence():
    """compute_steenrod_barcodes must not require a prior compute_persistence."""
    st = _rp2_simplex_tree()
    _, steenrod = st.compute_steenrod_barcodes(k=1)
    assert len(steenrod) > 2


def test_rejects_negative_k():
    with pytest.raises(Exception):
        _rp2_simplex_tree().compute_steenrod_barcodes(k=-1)


# ---------------------------------------------------------------------------
# CCP² — Sq², finite bar with integer filtration
#
# Filtration: CP² simplices ordered by (dim, lex) at values 0..254;
#             cone simplices ordered by (dim, lex) at values 255..510.
#
# N_CP2 = 255 is the boundary: indices 0..254 belong to CP²,
#                               indices 255..510 belong to the cone.
#
# With this filtration the bar lengths carry structural information:
#   ordinary[2]:  birth=73,  death=329, length=256  (H²(CP²))
#   ordinary[4]:  birth=254, death=510, length=256  (H⁴(CP²))
#   steenrod[4]:  birth=254, death=329, length=75   (Sq²: H²->H⁴)
#
# The Sq² bar is strictly shorter than the H⁴ bar (75 < 256): Sq² detects
# structure that ordinary persistence does not.
# ---------------------------------------------------------------------------

# 9-vertex Kühnel–Banchoff triangulation of CP² (36 top 4-simplices).
# Reference: Kühnel & Banchoff, Math. Intelligencer 5(3), 1983.
CP2_TOP = (
    (1, 2, 4, 5, 6), (2, 3, 5, 6, 4), (3, 1, 6, 4, 5),
    (1, 2, 4, 5, 9), (2, 3, 5, 6, 7), (3, 1, 6, 4, 8),
    (2, 3, 6, 4, 9), (3, 1, 4, 5, 7), (1, 2, 5, 6, 8),
    (3, 1, 5, 6, 9), (1, 2, 6, 4, 7), (2, 3, 4, 5, 8),
    (4, 5, 7, 8, 9), (5, 6, 8, 9, 7), (6, 4, 9, 7, 8),
    (4, 5, 7, 8, 3), (5, 6, 8, 9, 1), (6, 4, 9, 7, 2),
    (5, 6, 9, 7, 3), (6, 4, 7, 8, 1), (4, 5, 8, 9, 2),
    (6, 4, 8, 9, 3), (4, 5, 9, 7, 1), (5, 6, 7, 8, 2),
    (7, 8, 1, 2, 3), (8, 9, 2, 3, 1), (9, 7, 3, 1, 2),
    (7, 8, 1, 2, 6), (8, 9, 2, 3, 4), (9, 7, 3, 1, 5),
    (8, 9, 3, 1, 6), (9, 7, 1, 2, 4), (7, 8, 2, 3, 5),
    (9, 7, 2, 3, 6), (7, 8, 3, 1, 4), (8, 9, 1, 2, 5),
)
CONE_V = 10    # new vertex not in CP²
N_CP2  = 255   # number of faces in the CP² triangulation


def _cone_cp2_simplex_tree():
    """CCP² with integer filtration values.

    CP² faces sorted by (dim, lex) receive filtration values 0..N_CP2-1.
    Cone simplices — each CP² face joined with CONE_V, plus CONE_V alone
    (256 simplices total) — sorted by (dim, lex) receive values N_CP2..N_CP2+255.

    Note: cone simplices are built by joining CONE_V with each of the 255 CP²
    faces directly, not by taking all faces of the top cone simplices (which
    would incorrectly include CP² faces as sub-faces and corrupt the offsets).
    """
    cp2_faces = _all_faces(CP2_TOP)
    assert len(cp2_faces) == N_CP2

    cone_faces = {tuple(sorted(list(f) + [CONE_V])) for f in cp2_faces}
    cone_faces.add((CONE_V,))
    assert len(cone_faces) == 256

    st, _ = _simplex_tree_with_integer_filtration(cp2_faces, offset=0)
    for i, f in enumerate(sorted(cone_faces, key=lambda f: (len(f), f))):
        st.insert(list(f), filtration=float(N_CP2 + i))

    return st


def test_cone_cp2_total_simplex_count():
    """CCP² should have exactly 511 simplices (255 CP² + 256 cone)."""
    assert _cone_cp2_simplex_tree().num_simplices() == 511


def test_cone_cp2_sq2_one_finite_bar_in_dim4():
    """CCP² has exactly one Sq^2 bar in absolute cohomological degree 4.

    The Sq^2 map H^2(CP^2; F_2) -> H^4(CP^2; F_2) is an isomorphism, producing
    one finite persistence bar.  With integer filtration the bar is
    (birth=254, death=329), spanning 75 simplices.
    """
    _, steenrod = _cone_cp2_simplex_tree().compute_steenrod_barcodes(k=2, absolute=True)
    sq2 = steenrod[4]
    assert len(sq2) == 1
    birth, death = sq2[0]
    assert birth == 254.0   # last CP² 4-simplex completes H⁴(CP²)
    assert death == 329.0   # cone 3-simplex kills the Sq² class
    assert death - birth == 75.0


def test_cone_cp2_sq2_bar_born_in_cp2_killed_in_cone():
    """The Sq^2 bar must be born in the CP² part and killed in the cone part."""
    _, steenrod = _cone_cp2_simplex_tree().compute_steenrod_barcodes(k=2, absolute=True)
    birth, death = steenrod[4][0]
    assert birth < N_CP2   # born before the cone is attached
    assert death >= N_CP2  # killed by a cone simplex


def test_cone_cp2_sq2_shorter_than_h4_bar():
    """The Sq^2 bar (length 75) is strictly shorter than the H^4 bar (length 256).

    Ordinary persistence cannot distinguish CP^2 from S^2 v S^4 (both have the
    same Betti numbers); the Steenrod barcode can, precisely because Sq^2 on
    H^2(CP^2) is non-trivial while Sq^2 on H^2(S^2 v S^4) is zero.
    The shorter bar length is one manifestation of this richer structure.
    """
    ordinary, steenrod = _cone_cp2_simplex_tree().compute_steenrod_barcodes(k=2, absolute=True)
    # The topologically meaningful H^4 bar is the one that spans the CP²/cone
    # boundary (born in CP², killed by the cone).
    h4_boundary = [(b, d) for b, d in ordinary[4]
                   if not math.isinf(d) and b < N_CP2 and d >= N_CP2]
    sq2_bars = steenrod[4]
    assert len(h4_boundary) == 1
    assert len(sq2_bars) == 1
    h4_length  = h4_boundary[0][1] - h4_boundary[0][0]
    sq2_length = sq2_bars[0][1]    - sq2_bars[0][0]
    assert sq2_length < h4_length


def test_cone_cp2_ordinary_h0_essential_h2_h4_boundary_bar():
    """CCP² ordinary barcode structure.

    H^0: one essential bar (the single connected component survives).
    H^2 and H^4: the filtration produces many finite bars within CP² itself
    (from the fine integer ordering), but exactly one bar in each degree spans
    the CP²/cone boundary — born in the CP² part (index < N_CP2) and killed
    by a cone simplex (index >= N_CP2).  These are the topologically meaningful
    bars corresponding to H^2(CP^2) and H^4(CP^2).
    """
    ordinary, _ = _cone_cp2_simplex_tree().compute_steenrod_barcodes(k=2, absolute=True)

    h0_ess = [b for b in ordinary[0] if math.isinf(b[1])]
    assert len(h0_ess) == 1

    # Exactly one boundary-spanning bar in each of H^2 and H^4.
    for deg in (2, 4):
        boundary = [(b, d) for b, d in ordinary[deg]
                    if not math.isinf(d) and b < N_CP2 and d >= N_CP2]
        assert len(boundary) == 1, (
            f"Expected 1 boundary-spanning H^{deg} bar, got {len(boundary)}"
        )
        birth, death = boundary[0]
        assert birth < N_CP2
        assert death >= N_CP2


def test_cone_cp2_relative_convention():
    """With absolute=False the Sq^2 bar appears at index 5 with birth > death."""
    _, steenrod = _cone_cp2_simplex_tree().compute_steenrod_barcodes(k=2, absolute=False)
    sq2_rel = steenrod[5]
    assert len(sq2_rel) == 1
    birth, death = sq2_rel[0]
    assert birth > death   # relative convention: birth = killing simplex (later)
    assert birth == 329.0
    assert death == 254.0


def test_absolute_default_is_true():
    """absolute=True must be the default."""
    st = _cone_cp2_simplex_tree()
    _, steenrod_default  = st.compute_steenrod_barcodes(k=2)
    _, steenrod_explicit = st.compute_steenrod_barcodes(k=2, absolute=True)
    assert steenrod_default == steenrod_explicit
