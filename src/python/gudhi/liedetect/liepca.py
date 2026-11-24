"""
This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.

Author(s):       Henrique Ennes & Raphaël Tinarrage

Copyright (C) 2016 Inria & Institute of Science and Technology Austria
-----------------------------------------------------------------------------------------------------------------------

This module provides functions for dimension reduction and orthonormalization (Step 1 of LieDetect), as well as
Lie PCA (Step 2).

-----------------------------------------------------------------------------------------------------------------------

Dimension reduction:
    get_covariance_matrix
    print_covariance_eigenvalues
    project_on_minimal_subspace

Orthonormalization:
    print_norms
    batch_matrix_multiplication
    orthonormalize
    Orthonormalize

Lie PCA:
    get_lie_pca_operator

-----------------------------------------------------------------------------------------------------------------------
"""

# Standard imports.
from numbers import Integral, Real
from typing import Literal

# Third-party imports.
import numpy as np
import scipy

"""
-----------------------------------------------------------------------------------------------------------------------
Dimension reduction
-----------------------------------------------------------------------------------------------------------------------
"""


def get_covariance_matrix(
    pts: np.ndarray,
    center: bool = False,
    normalize: bool = False,
    orbit_dim: int | None = None,
) -> np.ndarray:
    """
    Computes the covariance matrix of a point cloud.

    Args:
        pts (np.ndarray): Array representing the points.
        center (bool, optional): If True, centers the points before computing the covariance. Defaults to False.
        normalize (bool, optional): If True, normalizes the covariance matrix by its Frobenius norm and rescales by
            sqrt(orbit_dim). This normalization makes it close to a projection matrix. Defaults to False.
        orbit_dim (int, optional): Dimension used for normalization. Required if normalize is True.

    Returns:
        np.ndarray: The covariance matrix of the points.
    """
    # Centers if needed
    if center:
        pts = pts - np.mean(pts, axis=0)

    # Computes covariance matrix
    cov = np.mean([np.outer(pt, pt) for pt in pts], axis=0)

    # Normalize if needed.
    if normalize:
        norm = np.linalg.norm(cov)
        if norm == 0 or not np.isfinite(norm):
            raise ValueError("Local covariance must have a finite, nonzero norm.")
        cov = cov / norm * np.sqrt(orbit_dim)

    return cov


def print_covariance_eigenvalues(pts) -> None:
    """
    Prints the eigenvalues of the covariance matrix of the given points
    in decreasing order.

    Args:
        pts(np.ndarray): Point cloud.
    """
    # Computes covariance matrix
    cov = get_covariance_matrix(pts)

    # Computes eigenvalues and normalize
    eigenvalues = np.sort(np.linalg.eigvals(cov).real)[::-1]
    eigenvalues = eigenvalues / np.sum(eigenvalues)

    print("Covariance eigenvalues:", *[f"{v:.1e} " for v in eigenvalues])


def project_on_minimal_subspace(
    pts: np.ndarray, threshold_eigenvalue: float
) -> np.ndarray:
    """
    Projects the points onto the minimal subspace they span, based on their covariance matrix. The projection
    is done by removing components with eigenvalues below a certain threshold, after normalization by L1 norm
    of eigenvalues.

    Args:
        pts (np.ndarray): Array of shape (nb_points, ambient_dim) representing
            the sampled points.
        threshold_eigenvalue (float, optional): Threshold for eigenvalues to
            consider a component significant.

    Returns:
        np.ndarray: Projected points in the minimal subspace.
    """
    # Computes covariance matrix
    cov = get_covariance_matrix(pts)

    # Gets eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(cov)

    # Normalizes eigenvalues by L1 norm
    eigenvalues = eigenvalues / np.sum(eigenvalues)

    # Filters out small eigenvalues
    significant_indices = np.where(eigenvalues > threshold_eigenvalue)[0]

    # Projects onto the minimal subspace
    projected_pts = pts @ eigenvectors[:, significant_indices]

    return projected_pts


"""
-----------------------------------------------------------------------------------------------------------------------
Orthonormalization
-----------------------------------------------------------------------------------------------------------------------
"""


def print_norms(pts: np.ndarray) -> None:
    """
    Prints the mean and standard deviation of the norms of the points. Useful to check normalization.

    Args:
        pts(np.ndarray): Point cloud.
    """
    norms = np.linalg.norm(pts, axis=1)
    mean_norm = np.mean(norms)
    std_norm = np.std(norms)

    print(f"Mean distance to origin: {mean_norm:.1e} ± {std_norm:.1e}")


def batch_matrix_multiplication(pts: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    """
    It will be used for simplifying some repeated operations down the road (DRY).

    Args:
        pts(np.ndarray): Batch of vectors to be multiplied by the matrix.
        matrix(np.ndarray): Matrix to multiply the vectors.

    Returns:
        (np.ndarray): The product of the matrix with each vector in the batch.
    """
    return pts @ matrix.T


def _as_real_points(pts: np.ndarray) -> np.ndarray:
    pts = np.asarray(pts)
    if pts.ndim != 2 or pts.shape[1] == 0 or pts.dtype.kind not in "iuf":
        raise ValueError("Points must be a two-dimensional array of real coordinates.")
    pts = np.asarray(pts, dtype=float)
    if not np.isfinite(pts).all():
        raise ValueError("Point coordinates must be finite.")
    return pts


def _whiten_centered(pts: np.ndarray, rank_tol: float | None):
    if rank_tol is None:
        if len(pts) <= pts.shape[1]:
            raise ValueError("Orthonormalization requires full-rank centered covariance.")
        cov = get_covariance_matrix(pts, center=True, normalize=False)
        eigenvalues, basis = scipy.linalg.eigh(cov)
        if eigenvalues[0] <= 0:
            raise ValueError("Orthonormalization requires full-rank centered covariance.")
        scales = np.sqrt(eigenvalues)
    else:
        if (
            isinstance(rank_tol, (bool, np.bool_))
            or not isinstance(rank_tol, Real)
            or not 0 <= rank_tol < 1
        ):
            raise ValueError("'rank_tol' must be a relative covariance threshold in [0, 1).")
        if len(pts) == 0:
            raise ValueError("Orthonormalization found no retained directions.")
        _, singular_values, right_vectors = scipy.linalg.svd(pts, full_matrices=False)
        cutoff = max(np.sqrt(rank_tol), np.finfo(pts.dtype).eps * max(pts.shape))
        if singular_values[0] == 0:
            raise ValueError("Orthonormalization found no retained directions.")
        retained = singular_values / singular_values[0] > cutoff
        if not retained.any():
            raise ValueError("Orthonormalization found no retained directions.")
        basis = right_vectors[retained].T
        scales = singular_values[retained] / np.sqrt(len(pts))

    if basis.shape[1] == pts.shape[1]:
        # Preserve the symmetric ambient-coordinate convention at full rank.
        forward = (basis / scales) @ basis.T
        inverse = (basis * scales) @ basis.T
    else:
        forward = (basis / scales).T
        inverse = basis * scales

    whitened = batch_matrix_multiplication(pts, forward)
    mean_norm = np.mean(np.linalg.norm(whitened, axis=1))
    return whitened / mean_norm, forward / mean_norm, inverse * mean_norm


def orthonormalize(
    pts: np.ndarray, *, rank_tol: float | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Center, whiten, and scale a finite real cloud to mean point norm one.

    With rank_tol=None, require full-rank centered covariance and keep
    ambient coordinates. Otherwise select covariance eigenvalues above
    rank_tol times the largest, subject to a numerical SVD floor.
    Return transformed points and the forward matrix, whose rows are
    the retained coordinates when reduction occurs.
    """
    pts = _as_real_points(pts)
    if len(pts) == 0:
        raise ValueError("Orthonormalization requires at least one point.")
    transformed, forward, _ = _whiten_centered(pts - np.mean(pts, axis=0), rank_tol)
    return transformed, forward


class Orthonormalize:
    """
    Orthonormalizes a point cloud by centering, homogenizing and normalizing it. This is Step 1 of LieDetect.

    Parameters:
        pts (np.ndarray): Finite real point cloud. Integer coordinates are
            converted to floating point.
        rank_tol (float, optional): Relative covariance eigenvalue threshold
            in [0, 1). Directions must exceed this fraction of the largest
            eigenvalue and the numerical SVD floor. The default None requires
            full-rank centered covariance.
    """

    def __init__(self, pts: np.ndarray, *, rank_tol: float | None = None) -> None:
        self.points = _as_real_points(pts).copy()
        self.rank_tol = rank_tol
        if len(self.points) == 0:
            raise ValueError("Orthonormalization requires at least one point.")
        self.mean = np.mean(self.points, axis=0)
        centered = self.points - self.mean
        self.cov = centered.T @ centered / len(centered)
        self.is_fitted = False

    def print_norms(self) -> None:
        """
        Prints the mean and standard deviation of the norms of the points. Useful to check normalization.
        """
        print_norms(self.points)

    def fit_transform(self) -> np.ndarray:
        """
        Fits the orthonormalization to the points and returns the orthonormalized points.

        Creates new attributes:
            - orth_points: Orthonormalized points.
            - orthonormal_transf: Matrix for the linear projection onto the orthonormal basis.

        Returns:
            pts_orth (np.ndarray): Orthonormalized points.
        """
        points = _as_real_points(self.points)
        if len(points) == 0:
            raise ValueError("Orthonormalization requires at least one point.")
        mean = np.mean(points, axis=0)
        centered = points - mean
        orth_pts, forward, inverse = _whiten_centered(centered, self.rank_tol)
        cov = centered.T @ centered / len(centered)
        self.mean = mean
        self.cov = cov
        self.orth_points = orth_pts
        self.orthonormal_transf = forward
        self._inverse_orthonormal_transf = inverse
        self.is_fitted = True

        return orth_pts

    def fit(self) -> None:
        """
        Only fits the orthonormalization to the points and creates new attributes:
            - orth_points: Orthonormalized points.
            - orthonormal_transf: The matrix for the linear
                projection onto orthonormalization.
        """
        _ = self.fit_transform()

    def transform(self, pts: np.ndarray) -> np.ndarray:
        """
        Applies the same orthonormalization transformation to a new set of points. Naturally assumes that either
        'fit' or 'fit_transform' has been previously called.

        Args:
            pts (np.ndarray): Point cloud.

        Returns:
            (np.ndarray): Whitened points in the fitted coordinate space.
        """
        if not self.is_fitted:
            raise RuntimeError(
                "'fit' or 'fit_transform' must be run before 'transform'."
            )

        pts = _as_real_points(pts)
        if pts.shape[1] != self.mean.size:
            raise ValueError(f"Points must have {self.mean.size} coordinates.")
        return batch_matrix_multiplication(pts - self.mean, self.orthonormal_transf)

    def inverse_transform(
        self,
        pts: np.ndarray | None = None,
        add_mean: bool = True,
        in_place: bool = False,
    ) -> np.ndarray:
        """
        Lift fitted coordinates back to the original ambient space.

        After dimension reduction, inverse_transform(transform(points))
        gives their orthogonal projection onto the retained affine subspace.
        Requires a previous call to fit or fit_transform.

        Args:
            pts (np.ndarray, optional): Point cloud to apply the inverse
                transform. If None, uses the transformed training points
                from the last fit. Defaults to None.
            add_mean (bool): If True, adds the training mean to the inverse
                transformation. Defaults to True.
            in_place (bool): If True, makes the object's points the inversed
                transformed points. Defaults to False.

        Returns:
            inversed_points (np.ndarray): Inversed transformed points.
        """
        if not self.is_fitted:
            raise RuntimeError("'fit' or 'fit_transform' must be run before 'inverse'.")

        if pts is None:
            pts = self.orth_points
        else:
            pts = _as_real_points(pts)
            fitted_dim = self.orthonormal_transf.shape[0]
            if pts.shape[1] != fitted_dim:
                raise ValueError(f"Points must have {fitted_dim} fitted coordinates.")

        inversed_points = batch_matrix_multiplication(pts, self._inverse_orthonormal_transf)
        if add_mean:
            inversed_points += self.mean

        if in_place:
            self.points = inversed_points

        return inversed_points


"""
-----------------------------------------------------------------------------------------------------------------------
Lie PCA
-----------------------------------------------------------------------------------------------------------------------
"""


def _compress_lie_pca_to_skew(lie_pca: np.ndarray, ambient_dim: int) -> np.ndarray:
    # Right and left multiplication by K applies skew projection to the
    # input and output, respectively. Flattening uses row-major order.
    indices = np.arange(ambient_dim**2)
    transposed = indices.reshape(ambient_dim, ambient_dim).T.ravel()
    right_projected = (lie_pca - lie_pca[:, transposed]) / 2
    compressed = (right_projected - right_projected[transposed, :]) / 2

    # Add S = (I + T) / 2, then remove roundoff asymmetry.
    compressed[indices, indices] += 0.5
    compressed[indices, transposed] += 0.5
    return (compressed + compressed.T) / 2


def get_lie_pca_operator(
    pts: np.ndarray,
    nb_neighbors: int,
    orbit_dim: int,
    method: Literal["PCA", "covariance"] = "PCA",
    compress_to_skew: bool = False,
    verbose: bool = False,
) -> np.ndarray:
    """
    Computes the Lie-PCA operator of a given point cloud. For the estimation of normal spaces, two options are
    available:
            - 'covariance': uses normalized local covariance matrices.
            - 'PCA': uses the leading right singular vectors of local
                differences. Each neighborhood must have numerical rank at
                least orbit_dim, using a relative cutoff of machine epsilon
                times the larger dimension of the difference matrix.
    The raw full-matrix operator is returned by default. With 'compress_to_skew',
    the operator is the identity on symmetric matrices and the compression of
    the raw operator on skew-symmetric matrices.

    Args:
        pts (np.ndarray): Array of shape (nb_points, ambient_dim) representing the sampled points on the orbit.
        nb_neighbors (int): Number of neighbors to use for tangent space estimation.
        orbit_dim (int): Dimension of the orbit for tangent space estimation.
        method (str, optional): Method to estimate normal spaces, either 'covariance' or 'PCA'. Defaults to 'PCA'.
        compress_to_skew (bool, optional): If True, return S + K L K for the raw
            operator L and the symmetric and skew projections S and K. Defaults to False.
        verbose (bool, optional): If True, prints eigenvalues and eigengap information. Defaults to False.

    Returns:
        np.ndarray: The Lie-PCA operator as a matrix of shape (ambient_dim**2, ambient_dim**2).
    """
    try:
        from sklearn.neighbors import KDTree
    except ModuleNotFoundError as exc:
        if exc.name != "sklearn":
            raise
        raise ImportError("Lie-PCA requires scikit-learn.") from exc

    pts = _as_real_points(pts)
    nb_points, ambient_dim = pts.shape
    if nb_points < 2:
        raise ValueError("Lie-PCA requires at least two points.")
    if (
        isinstance(nb_neighbors, (bool, np.bool_))
        or not isinstance(nb_neighbors, Integral)
        or not 1 <= nb_neighbors < nb_points
    ):
        raise ValueError("'nb_neighbors' must be an integer between 1 and the point count minus one.")
    if (
        isinstance(orbit_dim, (bool, np.bool_))
        or not isinstance(orbit_dim, Integral)
        or not 1 <= orbit_dim < ambient_dim
    ):
        raise ValueError("'orbit_dim' must be a positive integer strictly smaller than the ambient dimension.")
    if nb_neighbors < orbit_dim:
        raise ValueError("'nb_neighbors' must be at least 'orbit_dim'.")
    squared_norms = np.einsum("ij,ij->i", pts, pts)
    if not np.all(np.isfinite(squared_norms) & (squared_norms > 0)):
        raise ValueError("Lie-PCA points must have finite, nonzero squared norms.")

    if method not in ("PCA", "covariance"):
        raise ValueError(f"Method {method} not recognized.")

    kdt = KDTree(pts, leaf_size=nb_neighbors + 1, metric="euclidean")
    candidates = kdt.query(pts, nb_neighbors + 1, return_distance=False)
    identity = np.eye(ambient_dim)
    proj_normal_spaces = []
    for i, row in enumerate(candidates):
        neighbors = row[row != i][:nb_neighbors]
        differences = pts[i] - pts[neighbors]
        if method == "PCA":
            _, singular_values, right_vectors = scipy.linalg.svd(differences, full_matrices=False)
            if singular_values[0] == 0:
                raise ValueError(f"Local covariance must have a finite, nonzero norm at point {i}.")
            cutoff = np.finfo(differences.dtype).eps * max(differences.shape)
            numerical_rank = np.count_nonzero(singular_values / singular_values[0] > cutoff)
            if numerical_rank < orbit_dim:
                raise ValueError(
                    f"The PCA neighborhood at point {i} has numerical rank {numerical_rank}, "
                    f"below orbit_dim={orbit_dim}."
                )
            tangent_basis = right_vectors[:orbit_dim]
            tangent = tangent_basis.T @ tangent_basis
        else:
            tangent = get_covariance_matrix(
                differences, center=False, normalize=True, orbit_dim=orbit_dim
            )
        proj_normal_spaces.append(identity - tangent)

    # Compute projections on lines.
    proj_lines = [
        np.outer(pts[i, :], pts[i, :]) / squared_norms[i]
        for i in range(nb_points)
    ]
    # Entry ((a, b), (i, j)) is the mean of N[a, i] * P[j, b].
    size = ambient_dim**2
    lie_pca = np.einsum("nai,njb->abij", proj_normal_spaces, proj_lines).reshape(size, size)
    lie_pca /= nb_points

    if compress_to_skew:
        lie_pca = _compress_lie_pca_to_skew(lie_pca, ambient_dim)

    if verbose:
        vals = np.sort(np.linalg.eigvals(lie_pca).real)
        print("Lie PCA first eigenvalues:", *[f"{v:.1e} " for v in vals[:4]], end=" ")
        print(f"""\x1b[34mEigengap #{orbit_dim}:
            {(vals[orbit_dim] / vals[orbit_dim - 1]):.1e}\x1b[0m.""")

    return lie_pca
