"""
This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.

Author(s):       Henrique Ennes & Raphaël Tinarrage

Copyright (C) 2016 Inria & Institute of Science and Technology Austria
-----------------------------------------------------------------------------------------------------------------------

Implementation of the main class of the LieDetect module, OrbitFitter.

-----------------------------------------------------------------------------------------------------------------------
"""

# Standard imports.
from typing import Literal

# Third-party imports.
import numpy as np
import scipy

# Local imports.
from .liepca import get_lie_pca_operator
from .optimization import find_closest_algebra
from .utils import sample_from_lie_algebra


def print_hausdorff_distance(
    pts1: np.ndarray, pts2: np.ndarray, verbose: bool = True
) -> float:
    """
    Prints the Hausdorff distance from first point cloud to second.

    Args:
        pts1 (np.ndarray): First point cloud.
        pts2 (np.ndarray): Second point cloud.
        verbose (bool): If True, prints the distance. Defaults to True.

    Returns:
        hausdorff_dist (float): Estimated non-symmetric Hausdorff distance.
    """
    hausdorff_dist = scipy.spatial.distance.directed_hausdorff(pts1, pts2)[0]
    if verbose:
        print(f"""Non-symmetric \x1b[34mHausdorff distance:
              {hausdorff_dist:.3e}\x1b[0m.""")
    return hausdorff_dist


class OrbitFitter:
    """
    Basic class for LieDetect functionalities.

    Parameters:
        pts (np.ndarray): Point cloud of shape (n_points, ambient_dim).
    """

    def __init__(self, pts: np.ndarray) -> None:
        self.points = pts
        self.shape = pts.shape

        # Takes care of checking whether particular methods have been called.
        self.is_lie_pca = False
        self.is_closest_algebra = False
        self.is_sample_orbit = False

        # These will be updated as the methods get called.
        self.orbit_dim = None
        self.lie_pca_operator_ = None
        self.group = None
        self.representation_type_ = None
        self.algebra_ = None
        self.orbit_ = None
        self.hausdorff_distances_ = None

    """
    Step 2: LiePCA.
    """

    def lie_pca(
        self,
        nb_neighbors: int,
        orbit_dim: int,
        method: Literal["PCA", "covariance"] = "PCA",
        correction: bool = True,
        verbose: bool = False,
    ) -> np.ndarray:
        """
        Computes the LiePCA operator from the input point cloud.

        Args:
            nb_neighbors (int): Number of neighbors to consider for the Lie-PCA operator.
            orbit_dim (int): Dimension of the orbit to detect.
            method (str): Method to compute the Lie-PCA operator. Options are "PCA" (default) and "covariance".
            correction (bool): Whether to apply the bias correction to the Lie-PCA operator.
            verbose (bool): Whether to print progress and debug information.

        Returns:
            lie_pca (np.ndarray): The computed Lie-PCA operator.
        """
        self.orbit_dim = orbit_dim
        self.lie_pca_operator_ = get_lie_pca_operator(
            pts=self.points,
            nb_neighbors=nb_neighbors,
            orbit_dim=orbit_dim,
            method=method,
            correction=correction,
            verbose=verbose,
        )
        self.is_lie_pca = True

        return self.lie_pca_operator_

    def print_lie_pca_eigenvalues(self, return_vals: bool = False) -> np.ndarray | None:
        """
        Prints the eigenvalues of the Lie-PCA operator.

        Args:
            return_vals (bool): If True, returns the eigenvalues as
                a numpy array.

        Returns:
            None
        """
        # The following code goes a little against the principle of DRY,
        # but we will keep it for clarity
        # (and to avoid refracting too much :p).
        vals = np.sort(np.linalg.eigvals(self.lie_pca_operator_).real)
        print("Lie PCA first eigenvalues:", *[f"{v:.1e} " for v in vals[:4]], end=" ")

        if return_vals:
            return vals

    """
    Step 3: LieDetect.
    """

    def closest_algebra(
        self,
        group: str,
        group_dim: int | None = None,
        frequency_max: int | None = None,
        reps_to_test: list | None = None,
        span_ambient_space: bool = True,
        method: Literal["bottom_lie_pca", "full_lie_pca", "abelian"] = "bottom_lie_pca",
        verbose: bool = False,
        verbose_top_scores: bool = False,
    ) -> tuple[tuple, list[np.ndarray]]:
        """
        Finds the closest pushforward Lie algebra of a group (subspace of skew-symmetric matrices) to the given
        Lie PCA operator. This requires that `get_lie_pca` has been run beforehand.

        Args:
            group (str): The group ('torus', 'SU(2)', or 'SO(3)').
            group_dim (Optional[int]): Dimension of the torus (if group='torus'), otherwise ignored.
            frequency_max (int): Maximum frequency for lattice search (torus case).
            reps_to_test (list, optional): List of representations to test. If None, computed automatically.
            span_ambient_space (bool): Whether to restrict to representations with orbits spanning the ambient space.
                Defaults to True.
            method (str): Optimization method: 'bottom_lie_pca', 'full_lie_pca',  or 'abelian' (torus only). Defaults
                to 'bottom_lie_pca'.
            verbose (bool): If True, prints progress and results. Defaults to False.
            verbose_top_scores (bool): If True, prints top scoring representations. Defaults to False.

        Returns:
            optimal_rep: The optimal representation type found.
            optimal_algebra: List of matrices forming the optimal Lie algebra.
        """
        if not self.is_lie_pca:
            raise RuntimeError(
                "'get_lie_pca' must be run before finding the closest algebra."
            )

        self.group = group
        self.representation_type_, self.algebra_ = find_closest_algebra(
            group=self.group,
            lie_pca=self.lie_pca_operator_,
            group_dim=group_dim,
            frequency_max=frequency_max,
            reps_to_test=reps_to_test,
            span_ambient_space=span_ambient_space,
            method=method,
            verbose=verbose,
            verbose_top_scores=verbose_top_scores,
        )
        self.is_closest_algebra = True

        return self.representation_type_, self.algebra_

    """
    Step 4: Compute distance to orbit.
    """

    def sample_orbit(
        self,
        nb_points: int,
        method: Literal["evenly_spaced", "random_uniform"] = "evenly_spaced",
        verbose: bool = False,
        x: np.ndarray | None = None,
    ) -> np.ndarray:
        """
        Samples points on the orbit of a compact Lie group representation, given its Lie algebra generators. We suppose
        that the algebra is isomorphic to the canonical algebra estimate. This allows us to compute the periods, which
        are, otherwise, not stably computable from the algebra alone. In particular, it requires that
        `get_closest_algebra` has been run beforehand.

        Args:
            nb_points (int): Number of points to sample.
            method (str): Sampling method, 'evenly_spaced' or 'random_uniform'. Defaults to 'evenly_spaced'.
            verbose (bool): Whether to print information about the sampled orbit. Defaults to False.
            x (np.ndarray, optional): Initial vector to act on. If None, the first vector of the data is.
                Defaults to None.

        Returns:
            np.ndarray: Array of sampled points on the orbit.
        """
        if not self.is_closest_algebra:
            raise RuntimeError(
                "'get_closest_algebra' must be run before sampling an orbit."
            )

        if x is None:
            x = self.points[0]

        self.orbit_ = sample_from_lie_algebra(
            group=self.group,
            rep_type=self.representation_type_,
            algebra=self.algebra_,
            x=x,
            nb_points=nb_points,
            method=method,
            verbose=verbose,
        )
        self.is_sample_orbit = True

        self.hausdorff_distances_ = (
            print_hausdorff_distance(self.points, self.orbit_, verbose=False),
            print_hausdorff_distance(self.orbit_, self.points, verbose=False),
        )

        return self.orbit_
