"""
This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.

Author(s):       Henrique Ennes & Raphaël Tinarrage

Copyright (C) 2016 Inria & Institute of Science and Technology Austria
-----------------------------------------------------------------------------------------------------------------------


This module provides functions to sample points on orbits of representations, from their Lie algebra generators.
Namely, the Lie algebras are stored through bases (tuples of skew-symmetric matrices), that we suppose to be isomorphic
to the canonical pushforward algebras of the representations implemented in the module "algebra", through conjugation
by an orthogonal matrix. This assumption it crucial since, combined with the representation type of the algebra, we are
able to reconstruct the orbit. The evenly spaced torus grid uses joint periods of the connected
image group. Distinct transformations can give repeated points when the base point has a
nontrivial stabilizer in that image. Period calculations depend on the integer weight presentation.

-----------------------------------------------------------------------------------------------------------------------

Sample on orbits:
    sample_from_lie_algebra
    sample_from_lie_group_rep
    sample_from_lie_group
-----------------------------------------------------------------------------------------------------------------------
"""

# Standard imports.
from numbers import Integral
from typing import Literal

# Third-party imports.
import numpy as np
import scipy
from numpy.random import Generator, default_rng

RandomSeed = int | Generator | None

# Local imports.
from gudhi.liedetect.algebra import (
    get_canonical_pushforward_algebra,
    get_random_constrained_partition,
    get_random_lattice,
)
from gudhi.liedetect.utils import sample_from_lie_algebra

"""
-----------------------------------------------------------------------------------------------------------------------
Sample on orbits
-----------------------------------------------------------------------------------------------------------------------
"""


def _validate_base_point(x, ambient_dim):
    x = np.asarray(x)
    if x.shape != (ambient_dim,) or x.dtype.kind not in "iuf" or not np.isfinite(x).all():
        raise ValueError(f"'x' must be a finite real vector of length {ambient_dim}.")
    return x


def _is_cyclic(x, generators):
    # For a connected group, the orbit span is the generator-invariant
    # subspace obtained by repeatedly applying the generators to x.
    scale = np.max(np.abs(x.astype(float)))
    if scale == 0:
        return False
    initial = x / scale
    ambient_dim = x.size
    basis = np.empty((ambient_dim, ambient_dim))
    basis[:, 0] = initial / np.linalg.norm(initial)
    rank = 1
    tolerance = 10 * ambient_dim * np.finfo(basis.dtype).eps
    directions = [
        generator / norm
        for generator in generators
        if (norm := np.linalg.norm(generator)) > 0
    ]
    column = 0
    while column < rank and rank < ambient_dim:
        for direction in directions:
            residual = direction @ basis[:, column]
            accepted = basis[:, :rank]
            # Reorthogonalize before deciding whether this direction is new.
            for _ in range(2):
                residual -= accepted @ (accepted.T @ residual)
            norm = np.linalg.norm(residual)
            if norm > tolerance:
                basis[:, rank] = residual / norm
                rank += 1
                if rank == ambient_dim:
                    return True
        column += 1
    return rank == ambient_dim


def sample_from_lie_group_rep(
    group: str,
    rep_type: tuple,
    nb_points: int,
    x: np.ndarray | None = None,
    conjugate_algebra: bool = False,
    right_multiply_algebra: bool = False,
    translate_orbit: bool = False,
    method: Literal["evenly_spaced", "random_uniform"] = "evenly_spaced",
    verbose: bool = False,
    seed: RandomSeed = None,
) -> np.ndarray:
    """
    Samples points on the orbit of a compact Lie group representation, given its representation type. We suppose that
    the algebra is isomorphic to the canonical algebra indicated in rep_type. This allows us to compute the periods,
    which are, otherwise, not stably computable from the algebra alone.

    Args:
        group (str): The group type, e.g., 'torus' or 'SU(2)'.
        rep_type (tuple): Representation type parameters (e.g., weights or partition).
        nb_points (int): Number of points to sample.
        x (np.ndarray, optional): Initial vector to act on. Defaults to (1,1,...)/sqrt(ambient_dim).
        conjugate_algebra (bool, optional): Whether to conjugate the algebra by a random orthogonal matrix.
            Defaults to False.
        right_multiply_algebra (bool, optional): Whether to right-multiply the algebra by a random orthogonal matrix.
            Defaults to False.
        translate_orbit (bool, optional): Whether to rotate the sampled orbit by a random orthogonal transformation.
            Defaults to False.
        method (str, optional): Sampling method, 'evenly_spaced' or 'random_uniform'. Defaults to 'evenly_spaced'.
        verbose (bool, optional): Whether to print information about the sampled orbit.
        seed (int or numpy.random.Generator, optional): Seed or generator for sampling. If None, uses a random generator.

    Returns:
        np.ndarray: Array of sampled points on the orbit.
    """
    # Create generator.
    rng = default_rng(seed)
    # Validate the representation before reading its ambient dimension.
    algebra = get_canonical_pushforward_algebra(group=group, rep_type=rep_type)
    ambient_dim = algebra[0].shape[0]
    # Get initial vector.
    if x is None:
        x = np.ones(ambient_dim)
        x /= np.linalg.norm(x)
    else:
        x = _validate_base_point(x, ambient_dim)
    # Conjugate algebra if needed.
    if conjugate_algebra:
        orth = scipy.stats.special_ortho_group.rvs(ambient_dim, random_state=rng)
        algebra = [orth @ mat @ orth.T for mat in algebra]
        # Translate initial vector (to conserve a homogeneous orbit).
        x = orth @ x
    # Right-multiply algebra if needed.
    if right_multiply_algebra:
        if group == "torus":
            raise NotImplementedError(
                "Right-multiply not implemented for torus representations."
            )
        orth = scipy.stats.special_ortho_group.rvs(3, random_state=rng)
        algebra = [
            np.sum([algebra[j] * orth[j, i] for j in range(len(algebra))], axis=0)
            for i in range(len(algebra))
        ]
    # Sample orbit.
    orbit = sample_from_lie_algebra(
        group=group,
        rep_type=rep_type,
        algebra=algebra,
        nb_points=nb_points,
        x=x,
        method=method,
        verbose=verbose,
        seed=rng,
    )
    # Translate orbit if needed.
    if translate_orbit:
        orth = scipy.stats.special_ortho_group.rvs(ambient_dim, random_state=rng)
        orbit = np.array([orth @ point for point in orbit])
    return orbit


def sample_from_lie_group(
    group: str,
    ambient_dim: int,
    nb_points: int,
    frequency_max: int | None = None,
    group_dim: int | None = None,
    x: np.ndarray | None = None,
    conjugate_algebra: bool = False,
    right_multiply_algebra: bool = False,
    translate_orbit: bool = False,
    method: Literal["evenly_spaced", "random_uniform"] = "evenly_spaced",
    span_ambient_space: bool = False,
    verbose: bool = False,
    seed: RandomSeed = None,
) -> tuple[np.ndarray, tuple]:
    """
    Samples an orbit of a random representation in a given ambient space.

    Args:
        group (str): The group type, e.g., 'torus', 'SU(2)', or 'SO(3)'.
        ambient_dim (int): Dimension of the ambient space.
        nb_points (int): Number of points to sample on the orbit.
        frequency_max (Optional[int]): Maximal frequency for torus representations.
        group_dim (Optional[int]): Dimension of the group (for torus).
        x (np.ndarray, optional): Initial vector to act on. With non-Abelian full-span
            sampling and repeated irreps, the default is a seeded random unit vector.
            Otherwise the default is (1,1,...)/sqrt(ambient_dim).
        conjugate_algebra (bool, optional): Whether to conjugate the algebra by a random orthogonal matrix.
            Defaults to False.
        right_multiply_algebra (bool, optional): Whether to right-multiply the algebra by a random orthogonal matrix.
            Defaults to False.
        translate_orbit (bool, optional): Whether to rotate the sampled orbit by a random orthogonal transformation.
            Defaults to False.
        method (str, optional): Sampling method, 'random_uniform' or 'evenly_spaced'. Defaults to 'evenly_spaced'.
        span_ambient_space (bool, optional): Require the base point's continuous orbit
            to have full centered span. For tori, the nonzero character blocks must be
            distinct modulo sign and the base point must be nonzero in each plane.
            For SU(2) and SO(3), cyclicity is checked numerically and nearly dependent
            directions may be rejected. Defaults to False.
        verbose (bool, optional): Whether to print information about the sampled orbit. Defaults to False.
        seed (int or numpy.random.Generator, optional): Seed or generator for sampling. If None, uses a random generator.

    Returns:
        tuple[np.ndarray, tuple]: The sampled orbit (array of shape
            (nb_points, ambient_dim)) and representation type.
    """
    # Create generator.
    rng = default_rng(seed)
    # Gets random representation
    if group == "torus":
        if not isinstance(ambient_dim, Integral) or ambient_dim < 2 or ambient_dim % 2:
            raise ValueError("'ambient_dim' must be a positive even integer for a torus.")
        if span_ambient_space and x is not None:
            x = _validate_base_point(x, ambient_dim)
            if not np.all(np.any(x.reshape(-1, 2) != 0, axis=1)):
                raise ValueError(
                    "The torus base point must be nonzero in every rotation plane "
                    "when 'span_ambient_space' is True."
                )
        rep_type = get_random_lattice(
            lattice_rank=group_dim,
            ambient_rank=ambient_dim // 2,
            frequency_max=frequency_max,
            span_ambient_space=span_ambient_space,
            seed=rng,
        )

    elif group in ["SU(2)", "SO(3)"]:
        rep_type = get_random_constrained_partition(
            group=group,
            ambient_dim=ambient_dim,
            span_ambient_space=span_ambient_space,
            seed=rng,
        )
        if span_ambient_space:
            if x is None and len(set(rep_type)) != len(rep_type):
                x = rng.standard_normal(ambient_dim)
                x /= np.linalg.norm(x)
            if x is not None:
                x = _validate_base_point(x, ambient_dim)
                generators = get_canonical_pushforward_algebra(group, rep_type)
                if not _is_cyclic(x, generators):
                    raise ValueError(
                        f"The base point must be cyclic for the selected representation {rep_type} "
                        "when 'span_ambient_space' is True."
                    )

    else:
        raise NotImplementedError(f"Group '{group}' not recognized.")

    # Samples orbit from the representation type
    orbit = sample_from_lie_group_rep(
        group=group,
        rep_type=rep_type,
        nb_points=nb_points,
        x=x,
        conjugate_algebra=conjugate_algebra,
        right_multiply_algebra=right_multiply_algebra,
        translate_orbit=translate_orbit,
        method=method,
        verbose=verbose,
        seed=rng,
    )
    return orbit, rep_type
