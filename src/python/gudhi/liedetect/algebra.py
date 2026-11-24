"""
This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.

Author(s):       Henrique Ennes & Raphaël Tinarrage

Copyright (C) 2016 Inria & Institute of Science and Technology Austria
-----------------------------------------------------------------------------------------------------------------------

This module provides a collection of tools for linear algebra, lattice and partition computations, and canonical bases
of representations of the tori, SU(2), and SO(3). It includes utilities for manipulating skew-symmetric matrices.

-----------------------------------------------------------------------------------------------------------------------

Linear algebra:
    skew_sym_to_vect
    vect_to_skew_sym
    skew_sym_frame_to_projection
    gram_schmidt_orthonormalization

Lattices and partitions:
    get_random_lattice
    invariant_of_lattices
    get_lattices
    are_representations_equivalent
    get_partitions
    get_constrained_partitions
    get_random_constrained_partition

Canonical bases of representations:
    get_pushforward_alg_irrep_su2
    get_canonical_pushforward_algebra

-----------------------------------------------------------------------------------------------------------------------
"""

# Standard imports.
import sys
import itertools
from collections import Counter
from fractions import Fraction
from math import gcd
from numbers import Integral
from typing import Literal

# Third-party imports.
import numpy as np
from numpy.random import Generator, default_rng

RandomSeed = int | Generator | None

"""
-----------------------------------------------------------------------------------------------------------------------
Linear algebra
-----------------------------------------------------------------------------------------------------------------------
"""


def skew_sym_to_vect(mat: np.ndarray) -> np.ndarray:
    """
    Converts a skew-symmetric matrix mat, written as a matrix in the canonical basis of M_n(R), to its vector
    representation in the canonical basis of S_n(R). These are the matrices (-E_ij+E_ji) for i < n and i < j < n.
    In other words, the function returns the upper-diagonal entries of A as a vector.

    Args:
        mat (np.ndarray): A square skew-symmetric matrix of shape (n, n).

    Returns:
        np.ndarray: A vector of length n(n-1)/2 containing the sub-diagonal entries of A.

    Example:
        mat = np.array([[ 0,  1, 2],
                        [-1,  0, 1],
                        [-2, -1, 0]])
        vect = skew_sym_to_vect(mat)
        # vect = [1, 2, 1]
    """
    return _skew_sym_to_vect(mat, np)


def _skew_sym_to_vect(mat, array_backend):
    # Sanity check: the matrix must be skew-symmetric
    if not array_backend.allclose(mat, -mat.T):
        raise ValueError("Matrix A must be skew-symmetric.")

    # Gets the size of the matrix.
    n = np.shape(mat)[0]

    # Creates the indices of the canonical basis of S_n(R)
    indices = (tuple([i, j]) for i in range(n) for j in range(i + 1, n))

    # Gets entries at the indices
    vect = array_backend.array([mat[t] for t in indices])

    return vect


def vect_to_skew_sym(vect: np.ndarray) -> np.ndarray:
    """
    Converts a skew-symmetric matrix vect, written as a vector in the canonical basis of S_n(R), to its matrix
    representation in the canonical basis of M_n(R). In other word, the function reshapes the vector as an
    upper-diagonal matrix and skew-symmetrizes it.

    Args:
        vect (np.ndarray): mat vector of length n(n-1)/2 containing the sub-diagonal entries of a skew-symmetric
            matrix, ordered as (0,1), (0,2), ..., (0,n-1), (1,2), ..., (n-2,n-1).

    Returns:
        np.ndarray: The reconstructed skew-symmetric matrix of shape (n, n).

    Example:
        v = np.array([1, 2, 1])
        mat = vect_to_skew_sym(v)
        # mat =
        # [[ 0.  1.  2.]
        #  [-1.  0.  1.]
        #  [-2. -1.  0.]]
    """
    # Gets the size of the matrix
    n = int((np.sqrt(1 + 8 * len(vect)) + 1) / 2)

    # Creates the indices of the canonical basis of S_n(R)
    indices = ((i, j) for i in range(n) for j in range(i + 1, n))

    # Fills upper-diagonal entries of the matrix
    mat = np.zeros((n, n))
    for val, idx in zip(vect, indices):
        mat[idx] = val

    # Skew-symmetrize the matrix
    mat -= mat.T

    return mat


def skew_sym_frame_to_projection(
    frame: list[np.ndarray], method: Literal["QR", "differentiable"] = "QR"
) -> np.ndarray:
    """
    Given a list of d skew-symmetric matrices forming a frame (free family), returns the orthogonal projection matrix
    onto the subspace they span. This is an m x m matrix, where m is the dimension of S_n(R), i.e., m = n(n-1)/2.

    Args:
        frame (list[np.ndarray]): List of d skew-symmetric matrices of shape
            (n, n), assumed to be linearly independent.
        method (str, optional): Orthonormalization method.
            - "QR" uses NumPy's QR decomposition,
            - "differentiable" uses a manual Gram-Schmidt process
                (for autograd, where QR is not implemented).

    Returns:
        np.ndarray: The projection matrix of shape (m, m), where m = n(n-1)/2, representing the subspace in S_n(R).

    Example:
        mat = np.array([[ 0,  1, 2],
                      [-1,  0, 1],
                      [-2, -1, 0]])
        frame = [mat]
        projection = skew_sym_frame_to_projection(frame)
        # projection is a 3x3 projection matrix:
        # [[0.16666667 0.33333333 0.16666667]
        #  [0.33333333 0.66666667 0.33333333]
        #  [0.16666667 0.33333333 0.16666667]]
    """
    # Converts matrices to vectors in the canonical basis of S_n(R).
    if method == "differentiable":
        import autograd.numpy as array_backend
    else:
        array_backend = np
    frame_vectors = [_skew_sym_to_vect(matrix, array_backend) for matrix in frame]

    # Orthonormalizes
    frame_vectors = gram_schmidt_orthonormalization(frame_vectors, method=method)

    # Builds projection matrix.
    frame_vectors = array_backend.array(frame_vectors)
    projection = frame_vectors.T @ frame_vectors

    return projection


def gram_schmidt_orthonormalization(
    frame: list[np.ndarray], method: Literal["QR", "differentiable"] = "QR"
) -> list[np.ndarray]:
    """
    Orthonormalizes a list of vectors or square matrices via the Gram-Schmidt process.

    Args:
        frame (list[np.ndarray]): List of vectors or square matrices to
            orthonormalize.
        method (str, optional): Orthonormalization method.
            - "QR" uses NumPy's QR decomposition,
            - "differentiable" uses a manual Gram-Schmidt process
                (for autograd, where QR is not implemented).

    Returns:
        np.ndarray: List of orthonormalized vectors or matrices,
            matching the input type.
    """
    # Find whether the frame contains vectors or square matrices.
    is_matrix = frame[0].ndim > 1 and frame[0].shape[0] == frame[0].shape[1]
    n = frame[0].shape[0]

    if method == "differentiable":
        import autograd.numpy as array_backend
    else:
        array_backend = np

    # Normalize the frame.
    frame_orth = [v / array_backend.linalg.norm(v) for v in frame]

    # If matrix: flatten to vectors.
    if is_matrix:
        frame_orth = [mat.flatten() for mat in frame_orth]

    # Gram-Schmidt orthonormalization via np's QR decomposition
    if method == "QR":
        q, _ = np.linalg.qr(np.array(frame_orth).T)
        frame_orth = [v for v in q.T]

    # Gram-Schmidt orthonormalization process, manual implementation for autograd, where QR is not implemented
    if method == "differentiable":
        if len(frame_orth) > 1:
            for i in range(1, len(frame_orth)):
                v = frame_orth[i]
                for j in range(i):
                    w = frame_orth[j]
                    v -= array_backend.dot(v, w) * w
                v /= array_backend.linalg.norm(v)
                frame_orth[i] = v

    # If matrix: reshapes back to square matrices
    if is_matrix:
        frame_orth = [v.reshape((n, n)) for v in frame_orth]

    return frame_orth


"""
-----------------------------------------------------------------------------------------------------------------------
Lattices and partitions
-----------------------------------------------------------------------------------------------------------------------
"""


def get_random_lattice(
    lattice_rank: int,
    ambient_rank: int,
    frequency_max: int,
    span_ambient_space: bool = False,
    seed: RandomSeed = None,
) -> tuple[tuple[int, ...], ...]:
    """
    Generates a random lattice of rank lattice_rank in Z^ambient_rank. It may not be primitive.

    Args:
        lattice_rank (int): Dimension of the lattice.
        ambient_rank (int): Number of two-dimensional character blocks.
        frequency_max (int): Maximum frequency for irreps.
        span_ambient_space (bool): Whether to only consider representations whose orbits span the ambient space. Only
            implemented for rank-1 lattices.
        seed (int or numpy.random.Generator, optional): Seed or generator for sampling. If None, uses a random generator.

    Returns:
        tuple: The generated lattice as a tuple of tuples.
    """
    for name, value in (
        ("lattice_rank", lattice_rank),
        ("ambient_rank", ambient_rank),
        ("frequency_max", frequency_max),
    ):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 1:
            raise ValueError(f"'{name}' must be a positive integer.")
    if lattice_rank > ambient_rank:
        raise ValueError("Rank of ambient lattice is too small.")
    if span_ambient_space:
        if lattice_rank != 1:
            raise ValueError("'span_ambient_space' is only implemented for rank-1 lattices.")
        if ambient_rank > frequency_max:
            raise ValueError("Full span requires at least one distinct nonzero frequency per character block.")

    rng = default_rng(seed)
    for _ in range(100):
        if span_ambient_space:
            # Draw distinct nonzero magnitudes directly, then keep the existing primitive-row condition.
            magnitudes = rng.choice(frequency_max, size=ambient_rank, replace=False) + 1
            weights = tuple(int(value) for value in magnitudes * rng.choice((-1, 1), size=ambient_rank))
            if gcd(*weights) == 1:
                return (weights,)
        else:
            rows = rng.integers(-frequency_max, frequency_max, size=(lattice_rank, ambient_rank), endpoint=True)
            if np.linalg.matrix_rank(rows) == lattice_rank:
                return tuple(tuple(int(value) for value in row) for row in rows)
    raise RuntimeError("Failed to generate a valid lattice after 100 attempts.")


def invariant_of_lattices(
    lattice: tuple[tuple[int, ...], ...],
    method: Literal["span-equivalence", "orbit-equivalence"] = "span-equivalence",
    decimals_accuracy: int = 5,
) -> tuple:
    """
    Return the legacy numerical features for lattice scoring.

    These rounded features must not be used for exact representation equivalence.
        - 'span-equivalence':
            Returns the projection matrix on the space it spans.
        - 'orbit-equivalence':
            Returns the first projection matrix (for the lexicographic order) obtained by applying the Gram-Schmidt
            orthonormalization to all signed permutations.
    The argument "decimals_accuracy" is used for comparing the invariants of lattices.
    """
    # Special case: dimension 1.
    if len(lattice) == 1:
        return (tuple(np.sort(np.abs(lattice[0]))),)
    # General case for span-equivalence.
    elif method == "span-equivalence":
        frame = np.array(gram_schmidt_orthonormalization(np.asarray(lattice)))
        proj = np.sum([np.outer(v, v) for v in frame], axis=0)
        invariant = tuple(np.round(proj.flatten(), decimals=decimals_accuracy))
        return invariant
    # General case for orbit-equivalence.
    elif method == "orbit-equivalence":
        lattice = np.asarray(lattice)
        m = lattice.shape[1]
        return min(
            invariant_of_lattices(lattice[:, perm] * sign, "span-equivalence", decimals_accuracy)
            for perm in itertools.permutations(range(m))
            for sign in itertools.product([-1, 1], repeat=m)
        )
    else:
        raise ValueError(f"Method not recognized: {method}.")


def _rational_row_reduction(rows, pivot_columns):
    rows = [[Fraction(entry) for entry in row] for row in rows]
    rank = 0
    for column in range(pivot_columns):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        row = rows[rank]
        scale = row[column]
        row = rows[rank] = [entry / scale for entry in row]
        for i, other in enumerate(rows):
            if i != rank and other[column]:
                factor = other[column]
                rows[i] = [a - factor * b for a, b in zip(other, row)]
        rank += 1
        if rank == len(rows):
            break
    return tuple(tuple(row) for row in rows[:rank])


def _primitive_row(row):
    divisor = gcd(*row)
    if not divisor:
        return row
    if next(weight for weight in row if weight) < 0:
        divisor = -divisor
    return tuple(weight // divisor for weight in row)


def _exact_span(rows):
    if len(rows) == 1:
        return (_primitive_row(rows[0]),) if any(rows[0]) else ()
    return _rational_row_reduction(rows, len(rows[0]))


def _rational_projector(basis, width):
    rank = len(basis)
    if rank == 0 or rank == width:
        return tuple(tuple(int(rank != 0 and i == j) for j in range(width)) for i in range(width))
    # Solve (R R^T) X = R in exact arithmetic, then form R^T X.
    gram = [tuple(sum(a * b for a, b in zip(row, other)) for other in basis) for row in basis]
    solved = _rational_row_reduction(
        [gram_row + row for gram_row, row in zip(gram, basis)], rank
    )
    return tuple(
        tuple(sum(basis[k][i] * solved[k][rank + j] for k in range(rank)) for j in range(width))
        for i in range(width)
    )


def _projector_profiles(projector):
    return tuple(
        (row[i], tuple(sorted(map(abs, row)))) for i, row in enumerate(projector)
    )


def _projectors_orbit_equivalent(left, right, left_profiles, right_profiles):
    if left == right:
        return True
    width = len(left)
    choices = [
        [i for i, left_profile in enumerate(left_profiles) if left_profile == profile]
        for profile in right_profiles
    ]
    order = sorted(range(width), key=lambda index: len(choices[index]))
    permutation = [-1] * width
    used = set()

    def signs_are_consistent():
        # A nonzero edge prescribes the relative sign of its endpoints.
        # Propagate from one vertex per connected component and check cycles.
        signs = [0] * width
        for start in range(width):
            if signs[start]:
                continue
            signs[start] = 1
            pending = [start]
            while pending:
                i = pending.pop()
                for j in range(width):
                    if i == j or not right[i][j]:
                        continue
                    sign = signs[i] * (1 if right[i][j] == left[permutation[i]][permutation[j]] else -1)
                    if not signs[j]:
                        signs[j] = sign
                        pending.append(j)
                    elif signs[j] != sign:
                        return False
        return True

    def match(depth):
        if depth == width:
            return signs_are_consistent()
        j = order[depth]
        for i in choices[j]:
            if i in used or any(
                abs(left[i][permutation[k]]) != abs(right[j][k])
                for k in order[:depth]
            ):
                continue
            permutation[j] = i
            used.add(i)
            if match(depth + 1):
                return True
            used.remove(i)
        return False

    return match(0)


def get_lattices(
    lattice_rank: int,
    ambient_rank: int,
    frequency_max: int,
    method: Literal["span-equivalence", "orbit-equivalence"] = "span-equivalence",
    span_ambient_space: bool = True,
    verbose: bool = False,
) -> list[tuple[tuple[int, ...], ...]]:
    """
    Return full-rank integer row bases for torus candidate representations.

    With 'span-equivalence', retain the first representative of each rational
    row space in fixed coordinates. With 'orbit-equivalence', also identify
    row spaces related by signed column permutations. Representatives retain
    their original integer rows and can be nonprimitive.
    In addition, with the parameter "span_ambient_space", we only return lattices whose irreps in its decomposition do
    not repeat. This ensures that the corresponding orbit spans the ambient space.

    Args:
        ambient_rank (int): Dimension of the vector.
        lattice_rank (int): Rank of the lattice.
        frequency_max (int): Maximum frequency for irreps.
        method (str, optional): 'span-equivalence' or 'orbit-equivalence'.
        span_ambient_space (bool, optional): If True, only return lattices whose generic orbits span the ambient space.
        verbose (int, optional): Verbosity level.

    Returns:
        list: List of lattices, each as a tuple of frequency tuples.
    """
    if method not in ("span-equivalence", "orbit-equivalence"):
        raise ValueError(f"Method not recognized: {method}.")
    # Sanity check: ambient dimension must be large enough
    if not lattice_rank <= ambient_rank:
        raise ValueError("Rank of ambient lattice is too small.")

    if verbose:
        sys.stdout.write(f"""Generate lattices up to frequency {frequency_max}... """)

    # Gets irreps of T^lattice_rank with frequencies in
    # [0, frequency_max] (tuples of length lattice_dim)
    irreps = list(itertools.product(range(frequency_max + 1), repeat=lattice_rank))
    irreps.remove((0,) * lattice_rank)

    # Gets lattices, as combinations of irreps (comb instead of prod, for the orbit to span the ambient space)
    if span_ambient_space:
        lattices = itertools.combinations(irreps, ambient_rank)
    else:
        lattices = itertools.product(irreps, repeat=ambient_rank)

    # Visit candidates in the original order and keep the first representative of each span.
    full_rank_count = 0
    span_equivalence_classes = {}
    for characters in lattices:
        lattice = tuple(zip(*characters))
        if lattice_rank == 1 and gcd(*lattice[0]) != 1:
            continue
        span_key = _exact_span(lattice)
        if len(span_key) != lattice_rank:
            continue
        full_rank_count += 1
        span_equivalence_classes.setdefault(span_key, lattice)

    if verbose:
        sys.stdout.write(f"Full-rank lattices: {full_rank_count}... ")
        sys.stdout.write(f"Span-equivalence classes: {len(span_equivalence_classes)}. ")
    if method == "span-equivalence":
        return list(span_equivalence_classes.values())

    if lattice_rank == 1:
        rank_one_classes = {}
        for (row,), lattice in span_equivalence_classes.items():
            rank_one_classes.setdefault(tuple(sorted(map(abs, row))), lattice)
        representatives = list(rank_one_classes.values())
    else:
        representatives = []
        buckets = {}
        for basis, lattice in span_equivalence_classes.items():
            projector = _rational_projector(basis, ambient_rank)
            profiles = _projector_profiles(projector)
            bucket = buckets.setdefault(tuple(sorted(profiles)), [])
            if any(
                _projectors_orbit_equivalent(previous, projector, previous_profiles, profiles)
                for previous, previous_profiles in bucket
            ):
                continue
            bucket.append((projector, profiles))
            representatives.append(lattice)

    if verbose:
        sys.stdout.write(f"Orbit-equivalence classes: {len(representatives)}.\n")
    return representatives


def _torus_orbit_equivalent(rep0, rep1):
    _validate_torus_weights(rep0)
    _validate_torus_weights(rep1)
    width = len(rep0[0])
    if len(rep0) != len(rep1) or width != len(rep1[0]):
        return False
    left, right = (
        _exact_span(tuple(tuple(map(int, row)) for row in rep)) for rep in (rep0, rep1)
    )
    if len(left) != len(right):
        return False
    if left == right:
        return True
    if len(rep0) == 1:
        return sorted(map(abs, left[0])) == sorted(map(abs, right[0]))

    left_projector, right_projector = (
        _rational_projector(basis, width) for basis in (left, right)
    )
    left_profiles, right_profiles = (
        _projector_profiles(projector) for projector in (left_projector, right_projector)
    )
    return sorted(left_profiles) == sorted(right_profiles) and _projectors_orbit_equivalent(
        left_projector, right_projector, left_profiles, right_profiles
    )


def are_representations_equivalent(
    group: str, rep0: tuple, rep1: tuple, verbose: bool = False
) -> bool:
    """Check exact orbit-equivalence for a fixed source group.

    The connected image groups must be orthogonally conjugate. Torus inputs
    are nonempty rectangular integer weight tables. Their row counts specify
    the source torus rank and their column counts specify the ambient
    dimension. Dependent rows and zero columns are allowed. For SO(3) and
    SU(2), the complete irrep multisets, including trivial summands, must
    agree. With verbose=True, print the result.
    """
    if group in ("SU(2)", "SO(3)"):
        _validate_nonabelian_partition(group, rep0)
        _validate_nonabelian_partition(group, rep1)
        are_equivalent = sorted(rep0) == sorted(rep1)
    elif group == "torus":
        are_equivalent = _torus_orbit_equivalent(rep0, rep1)
    else:
        raise ValueError(f"Group not recognized: {group}.")

    if verbose:
        print(f"The representations {rep0} and {rep1}\x1b[1;31m are ", end="")
        if not are_equivalent:
            print("not ", end="")
        print("equivalent\x1b[0m.")
    return are_equivalent


def get_partitions(n: int):
    """Returns generator of partitions of n following
    https://jeromekelleher.net/generating-integer-partitions.html"""
    a = [0 for _ in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        ell = k + 1
        while x <= y:
            a[k] = x
            a[ell] = y
            yield a[: k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[: k + 1]


def get_constrained_partitions(
    group: str, ambient_dim: int, span_ambient_space: bool = False
) -> list[tuple[int, ...]]:
    """Return the nontrivial real representation partitions for the group.

    With span_ambient_space=True, each partition admits a cyclic vector
    and has no fixed directions, so it admits an orbit of full centered span.
    """
    # Defines integers (j % mod == rem) that are not irreps of the group
    if group == "SO(3)":
        mod, rem = 2, 0
    elif group == "SU(2)":
        mod, rem = 4, 2
    else:
        raise ValueError(f"Group {group} not recognized.")

    # Gets partitions satisfying the constraints
    partitions = []
    for partition in get_partitions(ambient_dim):
        fl = 1
        for j in partition:
            if j % mod == rem:
                fl = 0
                break
        if fl:
            partitions.append(tuple(partition))

    # Discards the trivial representation
    partitions.remove(tuple([1] * ambient_dim))

    if span_ambient_space:
        # Real-type irreps have commutant R. Even-dimensional SU(2) irreps
        # have commutant H, so their dimension over H is one quarter as large.
        partitions = [
            partition
            for partition in partitions
            if 1 not in partition
            and all(
                multiplicity <= (dim if dim % 2 else dim // 4)
                for dim, multiplicity in Counter(partition).items()
            )
        ]

    # Sanity check: the set cannot be empty
    if not partitions:
        raise ValueError("No partitions found.")
    return partitions


def get_random_constrained_partition(
    group: str,
    ambient_dim: int,
    span_ambient_space=False,
    seed: RandomSeed = None,
) -> tuple[int, ...]:
    """
    Returns a random partition of the integer n that is a valid representation of the specified group.
    """
    # Creates generator.
    rng = default_rng(seed)
    # Gets all partitions satisfying the constraints.
    partitions = get_constrained_partitions(
        group=group, ambient_dim=ambient_dim, span_ambient_space=span_ambient_space
    )
    # Select a random partition from the list.
    index = int(rng.integers(len(partitions)))
    return partitions[index]


"""
-----------------------------------------------------------------------------------------------------------------------
Canonical bases of representations
-----------------------------------------------------------------------------------------------------------------------
"""


def get_pushforward_alg_irrep_su2(dim: int) -> tuple[np.ndarray, ...] | None:
    """
    Returns a basis (x_1, x_2, x_3) of the pushforward Lie algebra for the irrep of SU(2) (and SO(3)) of dimension dim.
    """

    # Defines coefficients
    def delta(a, b):
        return (a == b) * 1

    def a_l(a, b):
        return np.sqrt((2 * a * b - a * (a - 1)) / 4)

    if dim == 1:
        j = 0
    elif dim % 2 == 1:
        j = int((dim - 1) / 2)
    elif dim % 4 == 0:
        j = (dim - 2) / 4
    else:
        print("Error:", dim, "is not a dimension of an irrep of SU(2).")
        return None

    # Defines matrices
    if type(j) is int:
        x_1 = np.zeros((int(2 * j) + 1, int(2 * j) + 1))
        x_2 = np.zeros((int(2 * j) + 1, int(2 * j) + 1))
        x_3 = np.zeros((int(2 * j) + 1, int(2 * j) + 1))
        for k in range(1, 2 * j + 2):
            for ell in range(1, 2 * j + 2):
                x_1[k - 1, ell - 1] = (
                    ((1 + (-1) ** k) / 2)
                    * (
                        delta(ell, k + 1) * a_l(int(k / 2), j)
                        + delta(ell + 3, k) * a_l(int((k - 2) / 2), j)
                    )
                    - (a_l(j, j) + np.sqrt((j**2 + j) / 2))
                    * (
                        delta(ell, 2 * j + 1) * delta(2 * j, k)
                        - delta(ell, 2 * j) * delta(2 * j + 1, k)
                    )
                    - ((1 + (-1) ** (k - 1)) / 2)
                    * (
                        delta(ell, k + 3) * a_l(int((k + 1) / 2), j)
                        + delta(ell + 1, k) * a_l(int((k - 1) / 2), j)
                    )
                )
                x_2[k - 1, ell - 1] = (
                    -(a_l(j, j) + np.sqrt((j**2 + j) / 2))
                    * (
                        delta(ell, 2 * j + 1) * delta(2 * j - 1, k)
                        - delta(ell, 2 * j - 1) * delta(2 * j + 1, k)
                    )
                    + delta(ell, k + 2) * a_l(int((k + 1) / 2), j)
                    - delta(ell + 2, k) * a_l(int((k - 1) / 2), j)
                )
                x_2[k - 1, ell - 1] = -x_2[k - 1, ell - 1]
                x_3[k - 1, ell - 1] = (
                    1
                    / 4
                    * (
                        (1 + (-1) ** k) * delta(ell + 1, k) * (2 * j + 2 - k)
                        + ((-1) ** k - 1) * delta(k + 1, ell) * (2 * j + 1 - k)
                    )
                )
    else:
        x_1 = np.zeros((int(4 * j) + 2, int(4 * j) + 2))
        x_2 = np.zeros((int(4 * j) + 2, int(4 * j) + 2))
        x_3 = np.zeros((int(4 * j) + 2, int(4 * j) + 2))
        for k in range(1, int(4 * j) + 3):
            for ell in range(1, int(4 * j) + 3):
                r = j
                x_1[k - 1, ell - 1] = ((1 + (-1) ** (k - 1)) / 2) * (
                    delta(ell, k + 3) * a_l(int((k + 1) / 2), r)
                    + delta(ell + 1, k) * a_l(int((k - 1) / 2), r)
                ) - ((1 + (-1) ** k) / 2) * (
                    delta(ell, k + 1) * a_l(int(k / 2), r)
                    + delta(ell + 3, k) * a_l(int((k - 2) / 2), r)
                )
                x_2[k - 1, ell - 1] = delta(ell, k + 2) * a_l(
                    int((k + 1) / 2), r
                ) - delta(ell + 2, k) * a_l(int((k - 1) / 2), r)
                x_3[k - 1, ell - 1] = (
                    1
                    / 4
                    * (
                        (1 + (-1) ** k) * delta(ell + 1, k) * (2 * j + 2 - k)
                        + ((-1) ** k - 1) * delta(k + 1, ell) * (2 * j + 1 - k)
                    )
                )
    return x_1, x_2, x_3


def _validate_torus_weights(rep_type):
    def is_row(row):
        return isinstance(row, (tuple, list, np.ndarray)) and np.ndim(row) == 1

    if (
        not isinstance(rep_type, (tuple, list, np.ndarray))
        or isinstance(rep_type, np.ndarray) and rep_type.ndim == 0
        or len(rep_type) == 0
        or not is_row(rep_type[0])
        or len(rep_type[0]) == 0
        or any(
            not is_row(row)
            or len(row) != len(rep_type[0])
            or any(isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) for value in row)
            for row in rep_type
        )
    ):
        raise ValueError("A torus representation requires a nonempty rectangular array of integer weights.")


def _validate_nonabelian_partition(group, rep_type):
    if (
        not isinstance(rep_type, (tuple, list, np.ndarray))
        or isinstance(rep_type, np.ndarray) and rep_type.ndim != 1
        or len(rep_type) == 0
        or any(
            isinstance(dim, (bool, np.bool_))
            or not isinstance(dim, Integral)
            or dim < 1
            or (dim % 2 == 0 and (group == "SO(3)" or dim % 4 != 0))
            for dim in rep_type
        )
    ):
        raise ValueError(f"Invalid real irrep dimensions for {group}: {rep_type!r}.")


def get_canonical_pushforward_algebra(group: str, rep_type: tuple) -> list[np.ndarray]:
    """
    Convert a representation-type of a Lie group into the canonical pushforward Lie algebra of the corresponding
    representation. The type can be
        - for the torus: a lattice basis,
        - for SU(2) and SO(3): a partition of an integer.
    """
    # If the group is the torus
    if group == "torus":
        _validate_torus_weights(rep_type)
        reduced_ambient_dim = len(rep_type[0])
        # Constructs basis of 2x2-block-diagonal skew-symmetric matrices
        basis = [
            np.zeros((2 * reduced_ambient_dim, 2 * reduced_ambient_dim))
            for _ in range(reduced_ambient_dim)
        ]
        for i in range(reduced_ambient_dim):
            basis[i][2 * i, 2 * i + 1], basis[i][2 * i + 1, 2 * i] = 1, -1
        # Generates infinitesimal rotations for the frequencies.
        pushforward_algebra = [
            np.sum([basis[j] * rep_type[i][j] for j in range(reduced_ambient_dim)], 0)
            for i in range(len(rep_type))
        ]

    # If the group is SU(2) or SO(3)
    elif group in ["SU(2)", "SO(3)"]:
        _validate_nonabelian_partition(group, rep_type)
        ambient_dim, nb_irreps = sum(rep_type), len(rep_type)
        algebra_irreps = [get_pushforward_alg_irrep_su2(k) for k in rep_type]
        pushforward_algebra = [np.zeros((ambient_dim, ambient_dim)) for _ in range(3)]
        index = 0
        for i in range(nb_irreps):
            k = rep_type[i]
            for j in range(3):
                pushforward_algebra[j][index : (index + k), index : (index + k)] = (
                    algebra_irreps[i][j]
                )
            index += k
    else:
        raise ValueError(f"Group not recognized: {group}.")
    return pushforward_algebra
