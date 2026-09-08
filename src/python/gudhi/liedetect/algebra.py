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
from math import gcd
from typing import Literal

# Third-party imports.
import autograd.numpy as np
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
    # Sanity check: the matrix must be skew-symmetric
    if not np.allclose(mat, -mat.T):
        raise ValueError("Matrix A must be skew-symmetric.")

    # Gets the size of the matrix.
    n = np.shape(mat)[0]

    # Creates the indices of the canonical basis of S_n(R)
    indices = (tuple([i, j]) for i in range(n) for j in range(i + 1, n))

    # Gets entries at the indices
    vect = np.array([mat[t] for t in indices])

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
    frame_vectors = [skew_sym_to_vect(matrix) for matrix in frame]

    # Orthonormalizes
    frame_vectors = gram_schmidt_orthonormalization(frame_vectors, method=method)

    # Builds projection matrix.
    frame_vectors = np.array(frame_vectors)
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

    # Normalize the frame.
    frame_orth = [v / np.linalg.norm(v) for v in frame]

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
                    v -= np.dot(v, w) * w
                v /= np.linalg.norm(v)
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
        ambient_rank (int): Ambient space dimension (should be even).
        frequency_max (int): Maximum frequency for irreps.
        span_ambient_space (bool): Whether to only consider representations whose orbits span the ambient space. Only
            implemented for rank-1 lattices.
        seed (int or numpy.random.Generator, optional): Seed or generator for sampling. If None, uses a random generator.

    Returns:
        tuple: The generated lattice as a tuple of tuples.
    """
    # Parameter checks.
    if frequency_max < 1:
        raise ValueError("'frequency_max' must be at least 1.")
    if lattice_rank > ambient_rank:
        raise ValueError("Rank of ambient lattice is too small.")
    # Create generator.
    rng = default_rng(seed)

    # Picks a random lattice
    has_maximal_rank = False

    max_attempts = 100
    for _ in range(max_attempts):
        # Generates lattice_rank random integral vectors in Z^ambient_rank
        values = list(range(-frequency_max, frequency_max + 1))
        lattice = tuple(
                    tuple(
                        int(value) for value in rng.choice(values, size=ambient_rank, replace=False, ))
                        for _ in range(lattice_rank)
        )
        # Check its rank
        has_maximal_rank = np.linalg.matrix_rank(np.array(lattice).T) == lattice_rank
        if has_maximal_rank:
            break
    else:
        raise RuntimeError("Failed to generate a valid lattice.")

    # If required, check whether the orbit spans the ambient space
    if span_ambient_space and lattice_rank == 1:
        # Checks whether the orbit spans the ambient space
        if (
            (0 in lattice[0])
            or (not gcd(*lattice[0]) == 1)
            or len(np.unique(np.abs(lattice[0]))) < len(lattice[0])
        ):
            # If not, generates a new lattice
            lattice = get_random_lattice(
                lattice_rank=lattice_rank,
                ambient_rank=ambient_rank,
                frequency_max=frequency_max,
                span_ambient_space=True,
                seed=rng,
            )

    elif span_ambient_space:
        raise ValueError("""The parameter 'span_ambient_space' is only
            implemented for rank-1 lattices.""")
    return lattice


def invariant_of_lattices(
    lattice: tuple[tuple[int, ...], ...],
    method: Literal["span-equivalence", "orbit-equivalence"] = "span-equivalence",
    decimals_accuracy: int = 5,
) -> tuple:
    """
    Returns an invariant for a lattice basis up to a certain equivalence relation.
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
        projections = []
        for perm in itertools.permutations(range(m)):
            for sign in itertools.product([-1, 1], repeat=m):
                permuted = lattice[:, perm] * sign  # apply signed permutation
                frame = np.array(gram_schmidt_orthonormalization(permuted))
                proj = np.sum([np.outer(v, v) for v in frame], axis=0)
                projections.append(np.round(proj.flatten(), decimals=decimals_accuracy))
        invariant = min(tuple(p) for p in projections)
        return invariant
    else:
        raise ValueError(f"Method not recognized: {method}.")


def get_lattices(
    lattice_rank: int,
    ambient_rank: int,
    frequency_max: int,
    method: Literal["span-equivalence", "orbit-equivalence"] = "span-equivalence",
    span_ambient_space: bool = True,
    verbose: bool = False,
) -> list[tuple[tuple[int, ...], ...]]:
    """
    Returns a list of lattices of rank lattice_dim in R^{ambient_rank}, written in a basis.
        - 'span-equivalence':
            Returns one lattice per span-equivalence class (span the same subspace of R^ambient_rank). This set is in
            correspondence with the primitive lattices, we do not guarantee that the selected representatives are
            primitive. For our purpose, having an arbitrary representative is enough.
        - 'orbit-equivalence':
            Returns one representative per orbit-equivalence class (up to signed permutations). It is stronger than
            'span-equivalence'.  For lattice_dim==1, the output are equivalent. Just as above, the output lattices may
            not be primitive.
    In addition, with the parameter "span_ambient_space", we only return lattices whose irreps in its decomposition do
    not repeat. This ensures that the corresponding orbit spans the ambient space.

    Args:
        ambient_rank (int): Dimension of the vector.
        lattice_rank (int): Rank of the lattice.
        frequency_max (int): Maximum frequency for irreps.
        method (str, optional): 'exact' for unique representatives, 'repetitions' for all (possibly repeated) lattices.
        span_ambient_space (bool, optional): If True, only return lattices whose generic orbits span the ambient space.
        verbose (int, optional): Verbosity level.

    Returns:
        list: List of lattices, each as a tuple of frequency tuples.
    """
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

    # Transposes the lattices, to write them as tuples of vectors
    # in Z^ambient_rank.
    lattices = [tuple(zip(*c)) for c in lattices]

    # For SO(2), discards non-primitive vectors
    if lattice_rank == 1:
        lattices = [lattice for lattice in lattices if gcd(*lattice[0]) == 1]

    # Keeps only maximal rank lattices
    lattices = [
        lattice
        for lattice in lattices
        if np.linalg.matrix_rank(lattice) == lattice_rank
    ]

    if verbose:
        sys.stdout.write(f"Full-rank lattices: {len(lattices)}... ")

    # Discards lattices that span the same vector subspace
    span_equivalence_classes = dict()
    for lattice in lattices:
        invariant = invariant_of_lattices(lattice=lattice, method="span-equivalence")
        if invariant not in span_equivalence_classes:
            span_equivalence_classes[invariant] = lattice
    lattices = list(span_equivalence_classes.values())

    if verbose:
        sys.stdout.write(f"Span-equivalence classes: {len(lattices)}"
                         f"{'... ' if method == 'orbit-equivalence' else '. '}")

    # If required, returns lattices obtained
    if method == "span-equivalence" or lattice_rank == 1:
        return lattices

    # Otherwise, discards orbit-equivalent lattices (under the action of signed permutations)
    elif method == "orbit-equivalence":
        orbit_equivalence_classes = {}
        for lattice in lattices:
            invariant = invariant_of_lattices(
                lattice=lattice, method="orbit-equivalence"
            )
            if invariant not in orbit_equivalence_classes:
                orbit_equivalence_classes[invariant] = lattice

        lattices = list(orbit_equivalence_classes.values())

        if verbose:
            sys.stdout.write(f"Orbit-equivalence classes: {len(lattices)}.\n")
        return lattices

    else:
        raise ValueError(f"Method not recognized: {method}.")


def are_representations_equivalent(
    group: str, rep0: tuple, rep1: tuple, verbose: bool = False
) -> bool:
    """
    Determines whether two representations are orbit-equivalent. For the non-Abelian groups SU(2) and SO(3), this boils
    down to the equality of irreps. For the torus, one checks whether the invariants of the lattices are equal.

    Args:
        group (str): The group type, one of "torus", "SU(2)", or "SO(3)". Defaults to "torus".
        rep0 (tuple): The first representation, as a tuple (e.g., partition or lattice basis).
        rep1 (tuple): The second representation, as a tuple.
        verbose (bool, optional): If True, prints the result. Defaults to False.

    Returns:
        bool: True if the representations are equivalent, False otherwise.
    """
    if group in ["SU(2)", "SO(3)"]:
        nonzero_irreps = [sorted([i for i in rep if i > 1]) for rep in [rep0, rep1]]
        are_equivalent = nonzero_irreps[0] == nonzero_irreps[1]
    elif group == "torus":
        invariants = [
            invariant_of_lattices(lattice=rep, method="orbit-equivalence")
            for rep in [rep0, rep1]
        ]
        are_equivalent = invariants[0] == invariants[1]
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
    """Returns all partitions of the integer n that are valid representations of the specified group."""
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

    # Discards representations that do not span the ambient space
    # (i.e., if contains the trivial irrep)
    if span_ambient_space:
        partitions = [
            partition
            for partition in partitions
            if (1 not in partition and len(np.unique(partition)) == len(partition))
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


def get_canonical_pushforward_algebra(group: str, rep_type: tuple) -> list[np.ndarray]:
    """
    Convert a representation-type of a Lie group into the canonical pushforward Lie algebra of the corresponding
    representation. The type can be
        - for the torus: a lattice basis,
        - for SU(2) and SO(3): a partition of an integer.
    """
    # If the group is the torus
    if group == "torus":
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
