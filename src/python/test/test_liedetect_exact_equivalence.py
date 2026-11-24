"""Independent exact orbit-equivalence checks for supported representations."""

import itertools
from math import gcd

import numpy as np
import pytest

from gudhi.liedetect import are_representations_equivalent
from gudhi.liedetect.algebra import get_lattices


def test_non_nested_coverings_give_the_same_orbit_class():
    covers = [((1, 2),), ((2, 4),), ((3, 6),)]
    for first, second in itertools.combinations(covers, 2):
        assert are_representations_equivalent("torus", first, second)
        assert are_representations_equivalent("torus", second, first)


def test_span_equivalence_preserves_the_weight_coordinate_order():
    assert get_lattices(1, 2, 2, "span-equivalence", False) == [
        ((1, 1),), ((1, 2),), ((2, 1),)
    ]
    assert get_lattices(1, 2, 2, "orbit-equivalence", False) == [
        ((1, 1),), ((1, 2),)
    ]


def test_close_integer_planes_are_not_merged():
    first = ((1, 0, 1000), (0, 1, 1000))
    second = ((1, 0, 1001), (0, 1, 1001))
    assert not are_representations_equivalent("torus", first, second)
    # Large integers can also alias if converted to float before elimination.
    large = 2**60
    assert not are_representations_equivalent(
        "torus", ((1, 0, large), (0, 1, large)), ((1, 0, large + 1), (0, 1, large + 1))
    )


def test_rational_rebasis_and_signed_permutations():
    rows = ((1, 0, 1, 0), (0, 1, 1, 2))
    # det([[2,1],[0,3]])=6, with a block permutation and a reflection.
    other = ((3, -2, 2, 1), (3, -6, 0, 3))
    assert are_representations_equivalent("torus", rows, other)


def test_source_rank_fixed_space_and_trivial_images():
    first = ((1, 2, 0), (2, 4, 0))
    second = ((0, 0, 0), (0, -2, 1))
    assert are_representations_equivalent("torus", first, second)
    assert not are_representations_equivalent("torus", first, ((1, 2, 0),))
    assert not are_representations_equivalent("torus", first, ((1, 2), (2, 4)))
    assert not are_representations_equivalent("torus", first, ((1, 2, 1), (2, 4, 2)))
    assert are_representations_equivalent("torus", ((0, 0),), ((0, 0),))
    assert not are_representations_equivalent("torus", ((0, 0),), ((1, 0),))


@pytest.mark.parametrize("group", ["SO(3)", "SU(2)"])
def test_nonabelian_partitions_include_fixed_directions(group):
    assert are_representations_equivalent(group, (1, 3, 3), (3, 1, 3))
    assert not are_representations_equivalent(group, (3,), (1, 3))
    assert not are_representations_equivalent(group, (1, 3, 3), (1, 1, 5))
    assert are_representations_equivalent(group, (1, 1), (1, 1))
    assert not are_representations_equivalent(group, (1,), (1, 1))


@pytest.mark.parametrize(
    ("group", "invalid"),
    [("torus", ((1.0, 2),)), ("torus", ((1,), (2, 3))), ("torus", ()), ("SO(3)", (4,))],
)
def test_invalid_representation_rejected(group, invalid):
    with pytest.raises(ValueError):
        are_representations_equivalent(group, invalid, invalid)


def _primitive_minors(rows):
    """Independent rank-two row-space key from Pluecker coordinates."""
    first, second = rows
    minors = tuple(
        first[i] * second[j] - first[j] * second[i]
        for i, j in itertools.combinations(range(len(first)), 2)
    )
    divisor = gcd(*minors)
    if next(value for value in minors if value) < 0:
        divisor = -divisor
    return tuple(value // divisor for value in minors)


def _exhaustive_orbit_key(rows):
    width = len(rows[0])
    return min(
        _primitive_minors(tuple(tuple(row[i] * sign for i, sign in zip(permutation, signs)) for row in rows))
        for permutation in itertools.permutations(range(width))
        for signs in itertools.product((-1, 1), repeat=width)
    )


def test_projector_matching_agrees_with_exhaustive_minor_oracle():
    rows = [
        ((1, 0, 1, 0), (0, 1, 0, 1)),
        ((1, 1, 0, 0), (0, 0, 1, -1)),
        ((1, 0, 1, 1), (0, 1, 1, -1)),
        ((1, 0, 1, 1), (0, 1, 1, 1)),
        ((1, 1, 1, 1), (1, -1, 0, 0)),
    ]
    signatures = [_exhaustive_orbit_key(rep) for rep in rows]
    for i, j in itertools.combinations(range(len(rows)), 2):
        assert are_representations_equivalent("torus", rows[i], rows[j]) == (
            signatures[i] == signatures[j]
        )


def test_equal_projector_diagonals_do_not_imply_orbit_equivalence():
    # Each table has orthogonal rows of equal norm, with all columns of
    # equal length. Both projectors have diagonal (1/2, 1/2, 1/2, 1/2).
    first = ((1, 0, 1, 0), (0, 1, 0, 1))
    second = ((5, 0, 3, -4), (0, 5, 4, 3))
    assert _exhaustive_orbit_key(first) != _exhaustive_orbit_key(second)
    assert not are_representations_equivalent("torus", first, second)
