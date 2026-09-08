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
                2 * pi / gcd(a1, ..., am) * 2 * ||(a1, ..., am)|| / ||A||.

    Case of T^d:
        Similar to the case of SO(2), but reasoning coordinate by coordinate.
    """

    def period_so2(weights: tuple[int, ...]) -> float:
        return 2 * np.pi / math.gcd(*weights)

    def norm_so2(weights: tuple[int, ...]) -> float:
        return np.sqrt(2) * np.linalg.norm(weights)

    # Compute the periods
    periods = [
        period_so2(weights) * norm_so2(weights) / np.linalg.norm(mat)
        for weights, mat in zip(rep_type, algebra)
    ]

    # Sanity check: the exponentiated matrices should be the identity
    for period, mat in zip(periods, algebra):
        if not np.isclose(scipy.linalg.expm(period * mat), np.eye(len(mat))).all():
            print(
                "Error! Incorrect period. Distance to identity:",
                np.linalg.norm(scipy.linalg.expm(period * mat) - np.eye(len(mat))),
            )
    return periods


def sample_orbit_from_algebra_torus(
    rep_type: tuple,
    algebra: list[np.ndarray],
    x: np.ndarray,
    nb_points: int,
    method: Literal["evenly_spaced", "random_uniform"] = "evenly_spaced",
    seed: RandomSeed = None,
) -> np.ndarray:
    """
    Sample points on the orbit of a torus representation from its Lie algebra. We suppose that the algebra is
    the canonical algebra indicated in rep_type.

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
    # Gets periods
    periods = np.asarray(get_periods_torus(rep_type, algebra), dtype=float)
    group_dim = len(algebra)

    # Generates grid based on the method.
    if method == "random_uniform":
        # Create generator.
        rng = default_rng(seed)
        # Generates random points in hypercube
        # [0, period[0]] x ... x [0, period[-1]]
        times = rng.random((nb_points, group_dim)) * periods
    elif method == "evenly_spaced":
        # Get number of points (potentially too many, sparsify later).
        nb_points_circle = int(np.ceil(nb_points ** (1.0 / group_dim)))
        grids = [
            np.linspace(0.0, periods[i], nb_points_circle, endpoint=False)
            for i in range(group_dim)
        ]
        times = np.array(list(itertools.product(*grids)), dtype=float)
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

    # Create generator.
    rng = default_rng(seed)

    # Create sampler.
    if method == "evenly_spaced":
        # With scramble=False, this is a fixed deterministic sequence.
        sampler = qmc.Halton(d=3, scramble=False)
        parameters = sampler.random(n=nb_points)
    elif method == "random_uniform":
        rng = np.random.default_rng(seed)
        parameters = rng.random((nb_points, 3))
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

    # Create orbit.
    orbit = np.empty((nb_points, x.size), dtype=float)
    for i in range(nb_points):
        group_element = (
            scipy.linalg.expm(alpha[i] * algebra_z)
            @ scipy.linalg.expm(beta[i] * algebra_y)
            @ scipy.linalg.expm(gamma[i] * algebra_z)
        )
        orbit[i] = group_element @ x
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
            rep_type=rep_type, algebra=algebra, x=x, nb_points=nb_points, seed=rng
        )
    else:
        raise ValueError(f"Group '{group}' not recognized.")
    if verbose:
        print(
            f"Sampled {len(orbit)} {method} points on the orbit of "
            f"\x1b[1;31m{group} with rep {rep_type}\x1b[0m."
        )
    return orbit
