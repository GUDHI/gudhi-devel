"""
This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.

Author(s):       Henrique Ennes & Raphaël Tinarrage

Copyright (C) 2016 Inria & Institute of Science and Technology Austria
-----------------------------------------------------------------------------------------------------------------------


This module provides util functions.

-----------------------------------------------------------------------------------------------------------------------

Sample on orbits:
    get_periods_torus
    sample_orbit_from_algebra_torus
    sample_orbit_from_algebra_su2
    sample_from_lie_algebra
----------------------------------------------------------------------------------------------------------------------
"""

# Standard imports.
import itertools
import math
from typing import Literal

# Third-party imports.
import numpy as np
import scipy
from numpy.random import Generator, default_rng
from scipy.stats import qmc

RandomSeed = int | Generator | None

# Local imports.
import gudhi.subsampling

from .algebra import _rational_row_reduction

"""
-----------------------------------------------------------------------------------------------------------------------
Sample on orbits
-----------------------------------------------------------------------------------------------------------------------
"""


def get_periods_torus(
    rep_type: tuple,
    algebra: list[np.ndarray],
) -> list[float]:
    """
    Returns a list of minimal periods t>0 such that exp(tA)=I for each A in the pushforward algebra. We suppose that
    the elements in "algebra" are such that they generate a periodic 1-parameter subgroup. More precisely, we suppose
    that the bijection between "algebra" and the canonical algebra (implemented in get_canonical_pushforward_algebra)
    is a Lie algebra isomorphism, induced by a conjugation. In particular, up to normalization by the norm of the
    elements, the periods are identical.

    Case of SO(2):
        The period of the integer matrix A* representing the pushforward algebra of SO(2) via the rep (a1, ..., am) is
                2 * pi / gcd(a1, ..., am).
        Its norm is
                sqrt(2) * ||(a1, ..., am)||.
        In particular, if A is a skew-symmetric matrix that is, up to normalization, conjugate to A*, then its period
            is
                2 * pi / gcd(a1, ..., am) * sqrt(2) * ||(a1, ..., am)|| / ||A||.

    Case of T^d:
        Similar to the case of SO(2), but reasoning coordinate by coordinate.

    A zero weight row acts as the identity at every time. It is assigned a
    sampling interval of length zero.
    """

    def period_so2(weights: tuple[int, ...]) -> float:
        return 2 * np.pi / math.gcd(*weights)

    def norm_so2(weights: tuple[int, ...]) -> float:
        return np.sqrt(2) * np.linalg.norm(weights)

    # Compute the periods
    periods = [
        period_so2(weights) * norm_so2(weights) / np.linalg.norm(mat) if any(weights) else 0.0
        for weights, mat in zip(rep_type, algebra)
    ]

    # Sanity check: the exponentiated matrices should be the identity
    for period, mat in zip(periods, algebra):
        if period and not np.isclose(scipy.linalg.expm(period * mat), np.eye(len(mat))).all():
            print(
                "Error! Incorrect period. Distance to identity:",
                np.linalg.norm(scipy.linalg.expm(period * mat) - np.eye(len(mat))),
            )
    return periods


def _integer_column_lattice_basis(weights):
    # Row Hermite elimination on W.T. Each operation is unimodular on
    # the generating rows, so the integer span is preserved.
    rows = [list(column) for column in zip(*weights)]
    rank = len(weights)
    for column in range(rank):
        pivot = next(i for i in range(column, len(rows)) if rows[i][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for i in range(column + 1, len(rows)):
            while rows[i][column]:
                quotient = rows[column][column] // rows[i][column]
                rows[column], rows[i] = rows[i], [
                    a - quotient * b for a, b in zip(rows[column], rows[i])
                ]
        if rows[column][column] < 0:
            rows[column] = [-value for value in rows[column]]
        for i in range(column):
            quotient = rows[i][column] // rows[column][column]
            rows[i] = [a - quotient * b for a, b in zip(rows[i], rows[column])]
    return tuple(zip(*rows[:rank]))


def _joint_period_generators(rep_type, algebra):
    # Retain original generator indices. Independent rows span the same
    # connected image even when the parameterization has redundant generators.
    weights = tuple(tuple(map(int, row)) for row in rep_type)
    independent = _rational_row_reduction(zip(*weights), len(weights))
    indices = [next(i for i, value in enumerate(row) if value) for row in independent]
    rank = len(indices)
    if rank == 0:
        return []
    independent_weights = [weights[i] for i in indices]

    basis = _integer_column_lattice_basis(independent_weights)
    augmented = [
        (*row, *(int(i == j) for j in range(rank)))
        for i, row in enumerate(basis)
    ]
    inverse = _rational_row_reduction(augmented, rank)
    # Convert the supplied generators back to the canonical angle scale.
    angle_generators = [
        algebra[i] * (
            np.sqrt(2) * np.linalg.norm(np.asarray(weights[i], dtype=float))
            / np.linalg.norm(algebra[i])
        )
        for i in indices
    ]
    # Row j of B^-1 supplies the coefficients of column j of B^-T.
    return [
        sum((2 * np.pi * float(coefficient)) * matrix
            for coefficient, matrix in zip(row[rank:], angle_generators))
        for row in inverse
    ]


def sample_orbit_from_algebra_torus(
    rep_type: tuple,
    algebra: list[np.ndarray],
    x: np.ndarray,
    nb_points: int,
    method: Literal["evenly_spaced", "random_uniform"] = "evenly_spaced",
    seed: RandomSeed = None,
) -> np.ndarray:
    """
    Sample a torus orbit using generators corresponding to the given integer rows.

    The generators may be simultaneously orthogonally conjugated and positively
    rescaled per row. The evenly spaced grid uses joint image periods. Distinct
    grid transformations give distinct orbit points when the base point has
    trivial stabilizer in the image group, up to floating-point precision.

    Args:
    rep_type (tuple): Representation type parameters (e.g., weights).
    algebra (list[np.ndarray]): List of Lie algebra generators as matrices.
    x (np.ndarray): Initial vector to act on.
    nb_points (int): Number of points to sample.
    method (str, optional): Sampling method, 'evenly_spaced' or 'random_uniform'. Defaults to 'evenly_spaced'.
    seed (int or numpy.random.Generator, optional): Seed or generator for sampling. If None, uses a random generator.

    Returns:
        np.ndarray: Array of sampled points on the orbit.
    """
    if method == "random_uniform":
        periods = np.asarray(get_periods_torus(rep_type, algebra), dtype=float)
        rng = default_rng(seed)
        times = rng.random((nb_points, len(algebra))) * periods
    elif method == "evenly_spaced":
        algebra = _joint_period_generators(rep_type, algebra)
        group_dim = len(algebra)
        if group_dim == 0:
            return np.repeat(np.asarray(x, dtype=float)[None, :], nb_points, axis=0)
        nb_points_circle = int(np.ceil(nb_points ** (1.0 / group_dim)))
        grid = np.linspace(0.0, 1.0, nb_points_circle, endpoint=False)
        times = np.array(list(itertools.product(grid, repeat=group_dim)), dtype=float)
    else:
        raise ValueError(
            f"Method '{method}' not recognized. Use 'evenly_spaced' or 'random_uniform'."
        )

    # Generates orbit (linear combinations or algebra elements wrt times)
    orbit = np.empty((len(times), x.size), dtype=float)
    for k, t_vec in enumerate(times):
        mat_alg = np.zeros_like(algebra[0], dtype=float)
        for t, mat in zip(t_vec, algebra):
            mat_alg += t * mat
        orbit[k] = scipy.linalg.expm(mat_alg) @ x
    # Sparsify if required.
    if len(orbit) > nb_points:
        orbit = gudhi.subsampling.choose_n_farthest_points(
            points=orbit,
            nb_points=nb_points,
            starting_point=0,
        )
    return np.array(orbit)


def sample_orbit_from_algebra_su2(
    rep_type: tuple,
    algebra: list[np.ndarray],
    x: np.ndarray,
    nb_points: int,
    method: Literal["evenly_spaced", "random_uniform"] = "evenly_spaced",
    seed: RandomSeed = None,
) -> np.ndarray:
    """
    Sample points on an orbit of an SU(2) or SO(3) representation.
    The group element is parameterized using the Euler decomposition:
            g(alpha, b, c) = exp(alpha Az) x exp(beta Ay) x exp(gamma Az)
    The "evenly_spaced" method uses a deterministic low-discrepancy sampling.
    The "random_uniform" method samples according to
        alpha ~ Uniform[0, 2π)
        beta ~ arccos(Uniform[-1, 1])
        gamma ~ Uniform[0, 4π) for SU(2) and Uniform[0, 2π) for SO(3)

    Args:
        rep_type (tuple): Representation type parameters (e.g., weights or partition).
        algebra (list[np.ndarray]): List of Lie algebra generators as matrices.
        x (np.ndarray): Initial vector to act on.
        nb_points (int): Number of points to sample.
        method (str, optional): Sampling method, 'evenly_spaced' or 'random_uniform'. Defaults to 'evenly_spaced'.
        seed (int or numpy.random.Generator, optional): Seed or generator for sampling. If None, uses a random generator.

    Returns:
        np.ndarray: Array of sampled points on the orbit.
    """

    def period_irrep_su2(dim: int) -> float:
        if dim % 4 == 0:
            return 4 * np.pi
        else:
            return 2 * np.pi

    # Gets algebra basis
    algebra_y = algebra[1]
    algebra_z = algebra[2]

    # Create sampler.
    if method == "evenly_spaced":
        # With scramble=False, this is a fixed deterministic sequence.
        sampler = qmc.Halton(d=3, scramble=False)
        parameters = sampler.random(n=nb_points)
    elif method == "random_uniform":
        parameters = default_rng(seed).random((nb_points, 3))
    else:
        raise ValueError(
            f"Method '{method}' not recognized. Use 'evenly_spaced' or 'random_uniform'."
        )

    # Sample coordinates.
    alpha = 2.0 * np.pi * parameters[:, 0]
    # Uniform cos(beta), not uniform beta.
    cos_beta = 2.0 * parameters[:, 1] - 1.0
    beta = np.arccos(np.clip(cos_beta, -1.0, 1.0))
    # For a representation factoring through SO(3), 2π suffices. SU(2) requires 4π.
    period = max(period_irrep_su2(dim) for dim in rep_type)
    gamma = period * parameters[:, 2]

    return _apply_euler_rotations(algebra_y, algebra_z, x, alpha, beta, gamma)


def _apply_euler_rotations(algebra_y, algebra_z, x, alpha, beta, gamma):
    # For real skew A, iA is Hermitian and exp(tA) is diagonal in its eigenbasis.
    eigenvalues_y, basis_y = scipy.linalg.eigh(1j * algebra_y)
    eigenvalues_z, basis_z = scipy.linalg.eigh(1j * algebra_z)
    coordinates = basis_z.conj().T @ x
    # Rows hold sample coordinates. This matrix changes them from the z basis to the y basis.
    z_to_y = basis_z.T @ basis_y.conj()
    y_to_z = z_to_y.conj().T

    orbit = np.empty((len(alpha), x.size), dtype=float)
    batch_size = 4096
    for start in range(0, len(orbit), batch_size):
        batch = slice(start, start + batch_size)
        vectors = np.exp(-1j * gamma[batch, None] * eigenvalues_z) * coordinates
        vectors = (vectors @ z_to_y) * np.exp(-1j * beta[batch, None] * eigenvalues_y)
        vectors = (vectors @ y_to_z) * np.exp(-1j * alpha[batch, None] * eigenvalues_z)
        orbit[batch] = (vectors @ basis_z.T).real
    return orbit


def sample_from_lie_algebra(
    group: str,
    rep_type: tuple,
    algebra: list[np.ndarray],
    nb_points: int,
    x: np.ndarray | None = None,
    method: Literal["evenly_spaced", "random_uniform"] = "evenly_spaced",
    verbose: bool = False,
    seed: RandomSeed = None,
) -> np.ndarray:
    """
    Samples points on the orbit of a compact Lie group representation, given its Lie algebra generators. We suppose
    that the algebra is isomorphic to the canonical algebra indicated in rep_type. This allows us to compute the
    periods, which are, otherwise, not stably computable from the algebra alone.

    Args:
        group (str): The group type, e.g., 'torus' or 'SU(2)'.
        rep_type (tuple): Representation type parameters (e.g., weights or partition).
        algebra (list[np.ndarray]): List of Lie algebra generators as matrices.
        nb_points (int): Number of points to sample.
        x (np.ndarray, optional): Initial vector to act on. Defaults to (1,1,...)/sqrt(ambient_dim).
        method (str, optional): Sampling method, 'evenly_spaced' or 'random_uniform'. Defaults to 'evenly_spaced'.
        verbose (bool, optional): Whether to print information about the sampled orbit.
        seed (int or numpy.random.Generator, optional): Seed or generator for sampling. If None, uses a random generator.

    Returns:
        np.ndarray: Array of sampled points on the orbit.
    """
    # Define base point.
    if x is None:
        ambient_dim = algebra[0].shape[0]
        x = np.ones(ambient_dim, dtype=float) / np.sqrt(ambient_dim)
    # Create generator.
    rng = default_rng(seed)
    # Generate orbit.
    if group == "torus":
        orbit = sample_orbit_from_algebra_torus(
            rep_type=rep_type,
            algebra=algebra,
            x=x,
            nb_points=nb_points,
            method=method,
            seed=rng,
        )
    elif group in ["SU(2)", "SO(3)"]:
        orbit = sample_orbit_from_algebra_su2(
            rep_type=rep_type,
            algebra=algebra,
            x=x,
            nb_points=nb_points,
            method=method,
            seed=rng,
        )
    else:
        raise ValueError(f"Group '{group}' not recognized.")
    if verbose:
        print(
            f"Sampled {len(orbit)} {method} points on the orbit of "
            f"\x1b[1;31m{group} with rep {rep_type}\x1b[0m."
        )
    return orbit
