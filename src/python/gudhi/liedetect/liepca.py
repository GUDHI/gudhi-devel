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
from typing import Literal

# Third-party imports.
import numpy as np
import scipy
import sklearn

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
        cov = cov / np.linalg.norm(cov) * np.sqrt(orbit_dim)

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
    return np.array([np.real(matrix.dot(x)) for x in pts])


def orthonormalize(pts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Orthonormalizes a point cloud by centering, homogenizing and normalizing it. This is Step 1 of LieDetect.

    Args:
        pts (np.ndarray): Point cloud.

    Returns:
        pts_orth (np.ndarray): Orthonormalized points.
        orthonormal_transf (np.ndarray): The matrix for the linear projection onto orthonormalization.
            Points are transformed by left multiplication by it.
    """
    # Copies to not modify the original points
    pts_orth = pts.copy()

    # Centers
    pts_orth -= np.mean(pts_orth, 0)

    # Homogenizes
    cov = get_covariance_matrix(pts_orth, center=True, normalize=False)
    orthonormal_transf = scipy.linalg.sqrtm(np.linalg.inv(cov))
    pts_orth = batch_matrix_multiplication(pts_orth, orthonormal_transf)

    # Normalizes
    mean_norm = np.mean(np.linalg.norm(pts_orth, axis=1))
    pts_orth /= mean_norm
    orthonormal_transf /= mean_norm

    return pts_orth, orthonormal_transf


class Orthonormalize:
    """
    Orthonormalizes a point cloud by centering, homogenizing and normalizing it. This is Step 1 of LieDetect.

    Parameters:
        pts (np.ndarray): Point cloud.
    """

    def __init__(self, pts: np.ndarray) -> None:
        self.points = pts.copy()
        self.cov = get_covariance_matrix(self.points, center=True, normalize=False)
        self.mean = np.mean(self.points, 0)
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
        orth_pts, orth_cov = orthonormalize(self.points)
        self.orth_points = orth_pts
        self.orthonormal_transf = orth_cov
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
            (np.ndarray): Point cloud projected by the same orthogonal
                transformation.
        """
        if not self.is_fitted:
            raise RuntimeError(
                "'fit' or 'fit_transform' must be run before 'transform'."
            )

        pts = pts.copy() - self.mean

        return batch_matrix_multiplication(pts, self.orthonormal_transf)

    def inverse_transform(
        self,
        pts: np.ndarray | None = None,
        add_mean: bool = True,
        in_place: bool = False,
    ) -> np.ndarray:
        """
        Inverses the orthonormalization. Naturally assumes that either 'fit' or 'fit_transform' has been previously
        called.

        Args:
            pts (np.ndarray, optional): Point cloud to apply the inverse
                transform. If None, applies to the training set. Defaults
                to None.
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
            pts = self.points

        inverse_transform = np.linalg.inv(self.orthonormal_transf)
        inversed_points = batch_matrix_multiplication(pts, inverse_transform)

        if in_place:
            self.points = inversed_points

        return inversed_points + add_mean * self.mean


"""
-----------------------------------------------------------------------------------------------------------------------
Lie PCA
-----------------------------------------------------------------------------------------------------------------------
"""


def get_lie_pca_operator(
    pts: np.ndarray,
    nb_neighbors: int,
    orbit_dim: int,
    method: Literal["PCA", "covariance"] = "PCA",
    correction: bool = True,
    verbose: bool = False,
) -> np.ndarray:
    """
    Computes the Lie-PCA operator of a given point cloud. For the estimation of normal spaces, two options are
    available:
            - 'covariance': uses local covariance matrices, correctly
                normalized.
            - 'PCA': uses local PCA, i.e., takes the top eigenvectors
                of previous estimation).
    In addition, the parameter 'correction' can be set to True to ensure that the operator is the identity on symmetric
    matrices (hence its kernel contains only skew-symmetric matrices).

    Args:
        pts (np.ndarray): Array of shape (nb_points, ambient_dim) representing the sampled points on the orbit.
        nb_neighbors (int): Number of neighbors to use for tangent space estimation.
        orbit_dim (int): Dimension of the orbit for tangent space estimation.
        method (str, optional): Method to estimate normal spaces, either 'covariance' or 'PCA'. Defaults to 'PCA'.
        correction (bool, optional): If True, applies a correction to ensure the operator is the identity on symmetric
            matrices. Defaults to True.
        verbose (bool, optional): If True, prints eigenvalues and eigengap information. Defaults to False.

    Returns:
        np.ndarray: The Lie-PCA operator as a matrix of shape (ambient_dim**2, ambient_dim**2).
    """
    nb_points, ambient_dim = np.shape(pts)

    # Computes local covariance matrices.
    kdt = sklearn.neighbors.KDTree(pts, leaf_size=nb_neighbors + 1, metric="euclidean")
    neighbors_idx = kdt.query(pts, nb_neighbors + 1, return_distance=False)[:, 1::]
    proj_tangent_spaces = [
        get_covariance_matrix(
            pts=pts[i] - pts[neighbors_idx[i, :]],
            center=False,
            normalize=True,
            orbit_dim=orbit_dim,
        )
        for i in range(nb_points)
    ]

    # Computes projection on normal spaces via local covariance (take complementary of previous estimation)
    if method == "covariance":
        proj_normal_spaces = [
            np.eye(ambient_dim) - proj for proj in proj_tangent_spaces
        ]

    # Computes projection on normal spaces via local PCA.
    elif method == "PCA":
        proj_normal_spaces = []
        for proj_tangent in proj_tangent_spaces:
            # Gets the eigenvalues and eigenvectors of the covariance matrix
            eigenvalues, mat = scipy.linalg.eigh(proj_tangent)

            # Gets top "orbit_dim" indices (largest eigenvalues), they represent the tangent space.
            idx = np.argsort(eigenvalues)[-orbit_dim:]

            # Creates canonical projection matrix (zero out the tangent space)
            proj_normal = np.eye(ambient_dim)
            proj_normal[idx, idx] = 0

            # Conjugates
            proj_normal_spaces.append(mat @ proj_normal @ mat.T)
    else:
        raise ValueError(f"Method {method} not recognized.")
    # Compute projections on lines.
    proj_lines = [
        np.outer(pts[i, :], pts[i, :]) / np.dot(pts[i, :], pts[i, :])
        for i in range(nb_points)
    ]
    # Create basis of space of matrices.
    basis_matrices = []
    for i in range(ambient_dim):
        for j in range(ambient_dim):
            mat = np.zeros((ambient_dim, ambient_dim))
            mat[i, j] = 1
            basis_matrices.append(mat)
    # Compute Lie-PCA operator.
    lie_pca = np.zeros((ambient_dim**2, ambient_dim**2))
    for i in range(len(basis_matrices)):
        lie_pca[:, i] = np.sum(
            [
                proj_normal_spaces[j] @ basis_matrices[i] @ proj_lines[j]
                for j in range(nb_points)
            ],
            axis=0,
        ).flatten()
    lie_pca /= len(pts)

    # Correction: set values of non-skew-symmetric matrices to zero. To do so, we skew-symmetrize the basis.
    if correction:
        lie_pca_corrected = np.zeros((ambient_dim**2, ambient_dim**2))
        for k in range(len(basis_matrices)):
            # Take basis element and decompose it into symmetric and skew-symmetric parts
            mat = basis_matrices[k]
            mat_sym = (mat + mat.T) / 2
            mat_skew_sym = (mat - mat.T) / 2

            # Computes the image via Lie PCA.
            # The image of the symmetric part is itself, so eigenvalue is 1
            mat_sym_image = mat_sym
            mat_skew_sym_image = (lie_pca @ (mat_skew_sym.reshape(-1))).reshape(
                ambient_dim, ambient_dim
            )
            mat_image = mat_sym_image + mat_skew_sym_image

            # Stores the image in the Lie PCA operator
            lie_pca_corrected[:, k] = mat_image.flatten()

        # Symmetrizes the operator
        lie_pca_corrected = (lie_pca_corrected + lie_pca_corrected.T) / 2
        lie_pca = lie_pca_corrected

    if verbose:
        vals = np.sort(np.linalg.eigvals(lie_pca).real)
        print("Lie PCA first eigenvalues:", *[f"{v:.1e} " for v in vals[:4]], end=" ")
        print(f"""\x1b[34mEigengap #{orbit_dim}:
            {(vals[orbit_dim] / vals[orbit_dim - 1]):.1e}\x1b[0m.""")

    return lie_pca
