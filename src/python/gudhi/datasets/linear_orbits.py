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
able to compute the periods of the basis elements (minimal t>0 such that exp(tA)=I), and hence be able to reconstruct
the orbit entirely (and without repetition for the torus). Note that, given the algebra alone, we cannot compute the
periods in a robust way.

-----------------------------------------------------------------------------------------------------------------------

Sample on orbits:
    sample_from_lie_algebra
    sample_from_lie_group_rep
    sample_from_lie_group
-----------------------------------------------------------------------------------------------------------------------
"""

# Standard imports.
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
        translate_orbit (bool, optional): Whether to translate the sampled orbit by a random orthogonal transformation.
            Defaults to False.
        method (str, optional): Sampling method, 'evenly_spaced' or 'random_uniform'. Defaults to 'evenly_spaced'.
        verbose (bool, optional): Whether to print information about the sampled orbit.
        seed (int or numpy.random.Generator, optional): Seed or generator for sampling. If None, uses a random generator.

    Returns:
        np.ndarray: Array of sampled points on the orbit.
    """
    # Create generator.
    rng = default_rng(seed)
    # Define ambient dimension.
    ambient_dim = len(rep_type[0]) * 2 if group == "torus" else sum(rep_type)
    # Get canonical pushforward algebra.
    algebra = get_canonical_pushforward_algebra(group=group, rep_type=rep_type)
    # Get initial vector.
    if x is None:
        x = np.ones(ambient_dim)
        x /= np.linalg.norm(x)
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
        x (np.ndarray, optional): Initial vector to act on. Defaults to (1,1,...)/sqrt(ambient_dim).
        conjugate_algebra (bool, optional): Whether to conjugate the algebra by a random orthogonal matrix.
            Defaults to False.
        right_multiply_algebra (bool, optional): Whether to right-multiply the algebra by a random orthogonal matrix.
            Defaults to False.
        translate_orbit (bool, optional): Whether to translate the sampled orbit by a random orthogonal transformation.
            Defaults to False.
        method (str, optional): Sampling method, 'random_uniform' or 'evenly_spaced'. Defaults to 'evenly_spaced'.
        span_ambient_space (bool, optional): Whether to only consider representations whose orbits span the ambient
            space. Defaults to False. Only implemented for the circle or the non-Abelian groups.
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
