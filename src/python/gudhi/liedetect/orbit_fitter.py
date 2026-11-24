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
        self.points = pts.copy()
        self.shape = self.points.shape
        self.is_lie_pca = False
        self.orbit_dim = None
        self.lie_pca_operator_ = None
        self._clear_algebra()

    def _clear_algebra(self) -> None:
        self.is_closest_algebra = False
        self.group = None
        self.representation_type_ = None
        self.algebra_ = None
        self.fit_diagnostics_ = None
        self._clear_orbit()

    def _clear_orbit(self) -> None:
        self.is_sample_orbit = False
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
        compress_to_skew: bool = False,
        verbose: bool = False,
    ) -> np.ndarray:
        """
        Compute the Lie-PCA operator used by all subsequent fitting methods.

        The default raw operator follows the full-matrix formulation in the
        LieDetect paper. Skew compression changes the finite-sample operator
        and the residuals used by the full fitting objective.

        Args:
            nb_neighbors (int): Number of neighbors to consider for the Lie-PCA operator.
            orbit_dim (int): Dimension of the orbit to detect.
            method (str): Method to compute the Lie-PCA operator. Options are "PCA" (default) and "covariance".
            compress_to_skew (bool): If True, use the skew compression of the raw
                operator and the identity on symmetric matrices. Defaults to False.
            verbose (bool): Whether to print progress and debug information.

        Returns:
            lie_pca (np.ndarray): The computed Lie-PCA operator.
        """
        operator = get_lie_pca_operator(
            pts=self.points,
            nb_neighbors=nb_neighbors,
            orbit_dim=orbit_dim,
            method=method,
            compress_to_skew=compress_to_skew,
            verbose=verbose,
        )
        self.orbit_dim = orbit_dim
        self.lie_pca_operator_ = operator
        self.is_lie_pca = True
        self._clear_algebra()

        return operator

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
        *,
        seed: int | np.random.Generator | None = None,
        initial_point: np.ndarray | None = None,
        optimizer_options: dict | None = None,
    ) -> tuple[tuple, list[np.ndarray]]:
        """
        Estimates the closest pushforward Lie algebra of a group (subspace of skew-symmetric matrices) to the given
        Lie PCA operator. This requires that `lie_pca` has been run beforehand.

        Args:
            group (str): The group ('torus', 'SU(2)', or 'SO(3)').
            group_dim (Optional[int]): Dimension of the torus (if group='torus'), otherwise ignored.
            frequency_max (int): Maximum frequency for lattice search (torus case).
            reps_to_test (list, optional): List of representations to test. If None, computed automatically.
            span_ambient_space (bool): Restrict the automatic search to representations
                admitting an orbit of full centered span. Defaults to True.
            method (str): Optimization method: 'bottom_lie_pca', 'full_lie_pca',  or 'abelian' (torus only). Defaults
                to 'bottom_lie_pca'.
            verbose (bool): If True, prints progress and results. Defaults to False.
            verbose_top_scores (bool): If True, prints top scoring representations. Defaults to False.
            seed (int or numpy.random.Generator, optional): Random source for optimizer starts.
            initial_point (np.ndarray, optional): Starting matrix in SO(ambient_dim) for each solver run.
                The second orthogonal component uses its reflected counterpart. Overrides the seed.
            optimizer_options (dict, optional): Pymanopt steepest descent settings overriding the defaults.
                These three controls apply to iterative solves.

        The successful fit sets `fit_diagnostics_` with the candidate objective value, stopping
        criterion, convergence flag, and per-run summaries. A direct solve has no iterative
        convergence flag. For `full_lie_pca`, the score is the sum of squared images
        of a Frobenius-orthonormal frame. A finite budget-limited result can participate
        in candidate selection.

        Returns:
            optimal_rep: The optimal representation type found.
            optimal_algebra: List of matrices forming the optimal Lie algebra.
        """
        if not self.is_lie_pca:
            raise RuntimeError(
                "'lie_pca' must be run before finding the closest algebra."
            )

        try:
            from .optimization import find_closest_algebra
        except ModuleNotFoundError as exc:
            if exc.name not in ("autograd", "pymanopt"):
                raise
            raise ImportError("LieDetect algebra fitting requires Autograd and Pymanopt.") from exc

        diagnostics = {}
        representation, algebra = find_closest_algebra(
            group=group,
            lie_pca=self.lie_pca_operator_,
            group_dim=group_dim,
            frequency_max=frequency_max,
            reps_to_test=reps_to_test,
            span_ambient_space=span_ambient_space,
            method=method,
            verbose=verbose,
            verbose_top_scores=verbose_top_scores,
            seed=seed,
            initial_point=initial_point,
            optimizer_options=optimizer_options,
            diagnostics=diagnostics,
        )
        self.group = group
        self.representation_type_ = representation
        self.algebra_ = algebra
        self.fit_diagnostics_ = diagnostics
        self.is_closest_algebra = True
        self._clear_orbit()

        return representation, algebra

    """
    Step 4: Compute distance to orbit.
    """

    def sample_orbit(
        self,
        nb_points: int,
        method: Literal["evenly_spaced", "random_uniform"] = "evenly_spaced",
        verbose: bool = False,
        x: np.ndarray | None = None,
        *,
        seed: int | np.random.Generator | None = None,
    ) -> np.ndarray:
        """
        Samples points on the orbit of a compact Lie group representation, given its Lie algebra generators. We suppose
        that the algebra is isomorphic to the canonical algebra estimate. This allows us to compute the periods, which
        are, otherwise, not stably computable from the algebra alone. In particular, it requires that
        `closest_algebra` has been run beforehand.

        Args:
            nb_points (int): Number of points to sample.
            method (str): Sampling method, 'evenly_spaced' or 'random_uniform'. Defaults to 'evenly_spaced'.
            verbose (bool): Whether to print information about the sampled orbit. Defaults to False.
            x (np.ndarray, optional): Initial vector to act on. If None, the first vector of the data is.
                Defaults to None.
            seed (int or numpy.random.Generator, optional): Random source for random sampling.

        Returns:
            np.ndarray: Array of sampled points on the orbit.
        """
        if not self.is_closest_algebra:
            raise RuntimeError(
                "'closest_algebra' must be run before sampling an orbit."
            )

        if x is None:
            x = self.points[0]

        orbit = sample_from_lie_algebra(
            group=self.group,
            rep_type=self.representation_type_,
            algebra=self.algebra_,
            x=x,
            nb_points=nb_points,
            method=method,
            verbose=verbose,
            seed=seed,
        )
        distances = (
            print_hausdorff_distance(self.points, orbit, verbose=False),
            print_hausdorff_distance(orbit, self.points, verbose=False),
        )
        self.orbit_ = orbit
        self.hausdorff_distances_ = distances
        self.is_sample_orbit = True

        return orbit
