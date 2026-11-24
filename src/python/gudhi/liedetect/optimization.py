"""
This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.

Author(s):       Henrique Ennes & Raphaël Tinarrage

Copyright (C) 2016 Inria & Institute of Science and Technology Austria
-----------------------------------------------------------------------------------------------------------------------

Implementation of the optimization methods to find the closest Lie algebra to a given Lie PCA operator, according to
three methods:
    "full_lie_pca": seeks a pushforward algebra that minimizes
        Lie PCA operator.
    "bottom_lie_pca": seeks a pushforward algebra close to the bottom
        eigenvectors of the Lie PCA operator.
    "abelian": computes the normal for of the bottom eigenvectors of the Lie
        PCA operator and finds the closest lattice.

-----------------------------------------------------------------------------------------------------------------------

Utils:
    chronometer_start
    chronometer_tick

Closest Lie algebra - general case:
    PARAMS_PYMANOPT
    find_closest_algebra
    optimization_bottom_lie_pca
    optimization_full_lie_pca

Closest Lie algebra - Abelian case:
    optimization_abelian
    normal_form_skew_sym_matrix
    normal_form_skew_symmetric_matrices

-----------------------------------------------------------------------------------------------------------------------
"""

# Standard imports
import datetime
import itertools
import sys
import time
from typing import Literal, Optional

# Third-party imports
import autograd.numpy as np
import pymanopt
import scipy
from numpy.random import default_rng
from pymanopt.optimizers.line_search import BackTrackingLineSearcher

# Local imports
from .algebra import (
    get_canonical_pushforward_algebra,
    get_constrained_partitions,
    get_lattices,
    invariant_of_lattices,
    skew_sym_frame_to_projection,
)

"""
-----------------------------------------------------------------------------------------------------------------------
Utils
-----------------------------------------------------------------------------------------------------------------------
"""


def chronometer_start(msg: str = "Start... ", verbose: bool = True) -> float:
    start_time = time.time()
    if verbose:
        sys.stdout.write(msg)
        sys.stdout.flush()
    return start_time


def chronometer_tick(start_time: float, i: int, i_total: int, msg: str) -> None:
    elapsed_time_secs = time.time() - start_time
    expected_time_secs = (i_total - i - 1) / (i + 1) * elapsed_time_secs
    msg1 = "It " + repr(i + 1) + "/" + repr(i_total) + ". "
    msg2 = "Duration %s. " % datetime.timedelta(seconds=round(elapsed_time_secs))
    msg3 = "Remaining %s." % datetime.timedelta(seconds=round(expected_time_secs))
    sys.stdout.write("\r" + msg + msg1 + msg2 + msg3)
    if i >= i_total - 1:
        sys.stdout.write("\n")


"""
-----------------------------------------------------------------------------------------------------------------------
Closest Lie algebra - general case
-----------------------------------------------------------------------------------------------------------------------
"""


PARAMS_PYMANOPT = {
    "max_iterations": 100,
    "min_gradient_norm": 1e-6,
    "max_time": 300,  # in seconds
    "verbosity": 0,
}


class _NonfiniteOptimizationError(RuntimeError):
    pass


class _StationaryBackTrackingLineSearcher(BackTrackingLineSearcher):
    def search(self, objective, manifold, x, d, f0, df0):
        if not (np.isfinite(f0) and np.isfinite(df0)):
            raise _NonfiniteOptimizationError("The cost or directional derivative is nonfinite.")
        # Steepest descent checks its stopping conditions after line search.
        if df0 == 0:
            return 0.0, x
        return super().search(objective, manifold, x, d, f0, df0)


def _run_optimizer(problem, ambient_dim, rng, initial_point, optimizer_options):
    options = PARAMS_PYMANOPT | (optimizer_options or {})
    if initial_point is None and rng is not None:
        initial_point = scipy.stats.special_ortho_group.rvs(
            ambient_dim,
            random_state=rng,
        )
    optimizer = pymanopt.optimizers.SteepestDescent(
        line_searcher=_StationaryBackTrackingLineSearcher(), **options
    )
    result = optimizer.run(problem, initial_point=initial_point)
    if not np.isfinite(result.cost) or not np.isfinite(result.point).all():
        raise _NonfiniteOptimizationError("The optimizer returned a nonfinite cost or point.")
    return result


def _run_summary(result, optimizer_options):
    tolerance = (optimizer_options or {}).get(
        "min_gradient_norm", PARAMS_PYMANOPT["min_gradient_norm"]
    )
    return {
        "cost": float(result.cost),
        "iterations": result.iterations,
        "stopping_criterion": result.stopping_criterion,
        "converged": bool(result.gradient_norm <= tolerance),
    }


def find_closest_algebra(
    group: str,
    lie_pca: np.ndarray,
    group_dim: Optional[int] = None,
    frequency_max: Optional[int] = None,
    reps_to_test=None,
    span_ambient_space: bool = True,
    method: Literal["bottom_lie_pca", "full_lie_pca", "abelian"] = "bottom_lie_pca",
    verbose: bool = True,
    verbose_top_scores: bool = False,
    *,
    seed=None,
    initial_point=None,
    optimizer_options=None,
    diagnostics=None,
) -> tuple[tuple, list[np.ndarray]]:
    """
    Finds the closest pushforward Lie algebra of a group (subspace of skew-symmetric matrices) to the given Lie PCA
    operator.

    Args:
        group (str): The group ('torus', 'SU(2)', or 'SO(3)').
        lie_pca (np.ndarray): The Lie PCA operator (matrix).
        group_dim (Optional[int]): Dimension of the torus (if group='torus'), otherwise ignored.
        frequency_max (int): Maximum frequency for lattice search (torus case).
        reps_to_test (list, optional): List of representations to test. If None, computed automatically.
        span_ambient_space (bool): Restrict the automatic search to representations admitting full centered span.
        method (str): Optimization method: 'bottom_lie_pca', 'full_lie_pca', or 'abelian' (torus only).
        verbose (bool): If True, print progress and results.
        verbose_top_scores (bool): If True, print top scoring representations.

    Returns:
        optimal_rep: The optimal representation type found.
        optimal_algebra: List of matrices forming the optimal Lie algebra.
    """
    if group not in ["torus", "SU(2)", "SO(3)"]:
        raise ValueError(f"Group {group!r} not recognized.")
    if method not in ("bottom_lie_pca", "full_lie_pca", "abelian"):
        raise ValueError(f"Method {method!r} not recognized.")
    if method == "abelian" and group != "torus":
        raise ValueError("The 'abelian' method is only supported for a torus.")
    group_dim = group_dim if group == "torus" else 3
    ambient_dim = int(np.sqrt(lie_pca.shape[0]))
    rng = None if seed is None else default_rng(seed)
    if initial_point is not None:
        initial_point = np.asarray(initial_point)
        if (
            initial_point.shape != (ambient_dim, ambient_dim)
            or initial_point.dtype.kind not in "iuf"
            or not np.isfinite(initial_point).all()
        ):
            raise ValueError("'initial_point' must be a finite matrix in SO(ambient_dim).")
        initial_point = initial_point.astype(float, copy=False)
        if (
            not np.allclose(initial_point.T @ initial_point, np.eye(ambient_dim))
            or np.linalg.det(initial_point) <= 0
        ):
            raise ValueError("'initial_point' must be a finite matrix in SO(ambient_dim).")
    runs = {}
    solver_kwargs = {
        "rng": rng,
        "initial_point": initial_point,
        "optimizer_options": optimizer_options,
    }


    if verbose:
        sys.stdout.write(f"""----> Optimization problem via \x1b[34m{method}\x1b[0m method for {group}  <----\n""")

    # Gets representations to test.
    explicit_candidates = reps_to_test is not None
    if reps_to_test is None:
        if group == "torus":
            if group_dim is None:
                raise ValueError("'group_dim' is required when group='torus'.")
            if frequency_max is None:
                raise ValueError("'frequency_max' is required when group='torus'.")
            if ambient_dim % 2 != 0:
                raise ValueError("'ambient_dim' must be even when group='torus'.")
            reps_to_test = get_lattices(
                lattice_rank=group_dim,
                ambient_rank=int(ambient_dim / 2),
                frequency_max=frequency_max,
                method="orbit-equivalence",
                span_ambient_space=span_ambient_space,
                verbose=verbose,
            )
        elif group in ["SU(2)", "SO(3)"]:
            reps_to_test = get_constrained_partitions(
                group=group,
                ambient_dim=ambient_dim,
                span_ambient_space=span_ambient_space,
            )
        else:
            raise ValueError(f"Group {group} not recognized.")

    if len(reps_to_test) == 0:
        raise ValueError("No candidate representations to test.")
    if explicit_candidates:
        if group == "torus" and group_dim is None:
            group_dim = len(reps_to_test[0])
        for rep_type in reps_to_test:
            algebra = get_canonical_pushforward_algebra(group, rep_type)
            if algebra[0].shape[0] != ambient_dim:
                raise ValueError(f"Candidate {rep_type!r} does not match ambient dimension {ambient_dim}.")
            if group == "torus" and (
                len(algebra) != group_dim or np.linalg.matrix_rank(rep_type) != group_dim
            ):
                raise ValueError(f"Torus candidate {rep_type!r} must have row rank {group_dim}.")

    # Gets bottom eigenvectors of Lie PCA if needed
    if method in ["bottom_lie_pca", "abelian"]:
        _, vecs = np.linalg.eigh(lie_pca)
        vecs = [
            vecs[:, i].reshape((ambient_dim, ambient_dim)) for i in range(group_dim)
        ]
        lie_pca_algebra = [(A - A.T) / 2 for A in vecs]
        lie_pca_proj = skew_sym_frame_to_projection(lie_pca_algebra)

    # Runs optimization in the general case
    if method in ["bottom_lie_pca", "full_lie_pca"]:
        if verbose:
            start_time = chronometer_start("Solve minimization problem... ")
        costs = dict()
        minimizers = dict()
        # Runs over all rep types
        for i in range(len(reps_to_test)):
            rep_type = reps_to_test[i]
            if method == "full_lie_pca":
                frame = _frobenius_orthonormal_frame(
                    get_canonical_pushforward_algebra(group, rep_type)
                )
            # Runs over the two connected components of the orthogonal group
            for determinant in ["+1", "-1"]:
                key = (rep_type, determinant)
                try:
                    if method == "bottom_lie_pca":
                        result = optimization_bottom_lie_pca(
                            group, lie_pca_proj, rep_type, determinant, **solver_kwargs
                        )
                    else:
                        result = _optimize_full_lie_pca_frame(
                            lie_pca, frame, determinant, **solver_kwargs
                        )
                except _NonfiniteOptimizationError as exc:
                    runs[key] = {"error": str(exc)}
                    continue
                costs[key] = result.cost
                minimizers[key] = result.point
                runs[key] = _run_summary(result, optimizer_options)
            if verbose:
                chronometer_tick(
                    start_time,
                    i,
                    len(reps_to_test),
                    "Solve minimization problem... ",
                )

        if not costs:
            raise RuntimeError("Algebra fitting failed: all optimizer runs were nonfinite.")
        winner = min(costs, key=costs.get)
        optimal_rep, _ = winner
        selected = runs[winner]
        optimal_change_of_basis = minimizers[winner]
        canonical_algebra = get_canonical_pushforward_algebra(
            group=group, rep_type=optimal_rep
        )
        optimal_algebra = [
            optimal_change_of_basis @ mat @ optimal_change_of_basis.T
            for mat in canonical_algebra
        ]

    # Runs optimization in the Abelian case
    elif method == "abelian":
        optimal_rep, optimal_algebra, costs = optimization_abelian(
            lie_pca_algebra=lie_pca_algebra, reps_to_test=reps_to_test,
            runs=runs, **solver_kwargs
        )
        winner = optimal_rep
        selected = runs.get(
            "joint_normal_form",
            {"iterations": 0, "stopping_criterion": "Direct Schur decomposition", "converged": None},
        )

    if diagnostics is not None:
        diagnostics.update(
            objective_value=float(costs[winner]),
            iterations=selected["iterations"],
            stopping_criterion=selected["stopping_criterion"],
            converged=selected["converged"],
            runs=runs,
        )

    if verbose:
        print(f"""The optimal rep found is \x1b[1;31m{optimal_rep}\x1b[0m with cost {min(costs.values()):.3e}.""")

    if verbose_top_scores:
        nb_scores_to_print = 10
        top_scores = sorted(costs.items(), key=lambda x: x[1])[:nb_scores_to_print]
        for i, (rep, score) in enumerate(top_scores):
            print(f"""    {rep} - cost {score:.3e}
                  (best cost #{i + 1}/{len(costs)})""")
    return optimal_rep, optimal_algebra


def optimization_bottom_lie_pca(
    group: str,
    lie_pca_proj: np.ndarray,
    rep_type: tuple,
    determinant: Literal["+1", "-1"] = "+1",
    *,
    rng=None,
    initial_point=None,
    optimizer_options=None,
):
    """
    Given a representation type and an initial guess of pushforward algebra, this function optimizes over the special
    orthogonal matrices to find a conjugation of the canonical matrices such that is the closest to the initial guess.
    We encode algebras (subspaces of the skew-symmetric matrices) as projection matrices.

    Args:
        group (str): The group type ('torus', 'SU(2)', or 'SO(3)').
        lie_pca_proj (np array): The target projection matrix on the space of skew-symmetric matrices.
        rep_type (tuple): The representation type
            (e.g., frequencies or partition).
        determinant (int): '+1' for SO(n), '-1' for the other component of O(n).

    Returns:
        result: A pymanopt optimization result with the optimal orthogonal
            matrix and cost.
    """
    # Defines basis
    canonical_algebra = get_canonical_pushforward_algebra(
        group=group, rep_type=rep_type
    )

    # Defines cost function.
    ambient_dim = len(canonical_algebra[0])
    manifold = pymanopt.manifolds.SpecialOrthogonalGroup(ambient_dim, k=1)

    @pymanopt.function.autograd(manifold)
    def _cost_function(orth: np.ndarray) -> float:
        # Transforms to provided connected component of orthogonal group
        if determinant == "-1":
            orth = orth @ np.diag([-1] + [1] * (ambient_dim - 1))
        # Conjugates canonical matrices
        algebra = [orth @ mat @ orth.T for mat in canonical_algebra]
        # Computes projection matrix
        proj = skew_sym_frame_to_projection(algebra, method="differentiable")
        # Computes distance to objective
        difference = proj - lie_pca_proj
        dist = np.trace(difference @ difference.T)
        return dist

    # Runs optimization
    problem = pymanopt.Problem(manifold, _cost_function)
    result = _run_optimizer(problem, ambient_dim, rng, initial_point, optimizer_options)
    # Transforms to provided connected component of orthogonal group
    if determinant == "-1":
        result.point = result.point @ np.diag([-1] + [1] * (ambient_dim - 1))
    return result


def _frobenius_orthonormal_frame(generators: list[np.ndarray]) -> np.ndarray:
    columns = np.stack([generator.ravel() for generator in generators], axis=1)
    basis, singular_values, _ = np.linalg.svd(columns, full_matrices=False)
    tolerance = np.finfo(columns.dtype).eps * max(columns.shape) * singular_values[0]
    if len(singular_values) != len(generators) or singular_values[-1] <= tolerance:
        raise ValueError("Full Lie-PCA requires linearly independent candidate generators.")
    return basis.T.reshape(len(generators), *generators[0].shape)


def optimization_full_lie_pca(
    group: str,
    lie_pca: np.ndarray,
    rep_type: tuple,
    determinant: Literal["+1", "-1"] = "+1",
    *,
    rng=None,
    initial_point=None,
    optimizer_options=None,
):
    """Fit a conjugation using the full Lie-PCA operator.

    The cost is the sum of squared images of a Frobenius-orthonormal frame
    under the operator in flattened matrix coordinates. The result contains
    the minimizing change of basis in the requested component of O(n).
    """
    frame = _frobenius_orthonormal_frame(get_canonical_pushforward_algebra(group, rep_type))
    return _optimize_full_lie_pca_frame(
        lie_pca,
        frame,
        determinant,
        rng=rng,
        initial_point=initial_point,
        optimizer_options=optimizer_options,
    )


def _optimize_full_lie_pca_frame(
    lie_pca, frame, determinant, *, rng, initial_point, optimizer_options
):
    ambient_dim = frame.shape[1]
    manifold = pymanopt.manifolds.SpecialOrthogonalGroup(ambient_dim, k=1)
    reflection = np.diag([-1] + [1] * (ambient_dim - 1))

    @pymanopt.function.autograd(manifold)
    def _cost_function(orth):
        if determinant == "-1":
            orth = orth @ reflection
        residuals = [lie_pca.dot((orth @ matrix @ orth.T).ravel()) for matrix in frame]
        return sum(residual @ residual for residual in residuals)

    problem = pymanopt.Problem(manifold, _cost_function)
    result = _run_optimizer(problem, ambient_dim, rng, initial_point, optimizer_options)
    if determinant == "-1":
        result.point = result.point @ reflection
    return result


"""
-----------------------------------------------------------------------------------------------------------------------
Closest Lie algebra - Abelian case
-----------------------------------------------------------------------------------------------------------------------
"""


def optimization_abelian(
    lie_pca_algebra: list, reps_to_test: list, *, rng=None,
    initial_point=None, optimizer_options=None, runs=None,
) -> tuple[tuple, list, dict]:
    """
    Finds the closest Abelian (torus) Lie algebra to the given Lie PCA algebra.

    Args:
        lie_pca_algebra (list): List of skew-symmetric matrices representing the Lie PCA algebra.
        reps_to_test (list): List of candidate lattice representations (as arrays) to compare.

    Returns:
        optimal_rep: The optimal lattice representation found (as a tuple of tuples).
        optimal_algebra: List of matrices forming the optimal Lie algebra in the original basis.
        costs: Dictionary mapping each tested lattice (as tuple of tuples) to its cost.
    """
    group_dim = len(lie_pca_algebra)
    # Finds optimal weights (as real numbers) via an optimization
    if group_dim == 1:
        opt_weights, opt_change_of_basis = normal_form_skew_sym_matrix(
            lie_pca_algebra[0]
        )
        opt_weights = (opt_weights,)
    else:
        opt_weights, opt_change_of_basis = normal_form_skew_symmetric_matrices(
            lie_pca_algebra, rng=rng, initial_point=initial_point,
            optimizer_options=optimizer_options, runs=runs,
        )

    # Computes invariant for the given lattice
    optimal_invariant = invariant_of_lattices(
        opt_weights,
        method="span-equivalence",
        decimals_accuracy=20,
    )
    # Compares to invariant of all the lattices to test
    # To allow the comparison, we take each lattice (it represents an
    # orbit-equivalence class) and compute its invariant under all
    # permutations of the columns and signs
    m = len(opt_weights[0])
    costs = dict()
    for weights in reps_to_test:
        if group_dim == 1:
            invariant = invariant_of_lattices(
                weights,
                method="span-equivalence",
                decimals_accuracy=20,
            )
            invariant /= np.linalg.norm(invariant)
            cost = np.linalg.norm(np.array(optimal_invariant) - np.array(invariant))
            if np.isfinite(cost):
                costs[weights] = cost
        else:
            lattice = np.asarray(weights)
            for perm in itertools.permutations(range(m)):
                for sign in itertools.product([-1, 1], repeat=m):
                    permuted = lattice[:, perm] * sign
                    invariant = invariant_of_lattices(
                        permuted,
                        method="span-equivalence",
                        decimals_accuracy=20,
                    )
                    permuted = tuple(tuple(row) for row in permuted.tolist())
                    cost = np.linalg.norm(np.array(optimal_invariant) - np.array(invariant))
                    if np.isfinite(cost):
                        costs[permuted] = cost

    # Deduces optimal representation
    if not costs:
        raise RuntimeError("Abelian algebra fitting failed: no finite candidate scores.")
    optimal_rep = min(costs, key=costs.get)

    # Defines optimal_algebra
    optimal_algebra = get_canonical_pushforward_algebra(
        group="torus", rep_type=optimal_rep
    )
    optimal_algebra = [
        opt_change_of_basis @ mat @ opt_change_of_basis.T for mat in optimal_algebra
    ]

    return optimal_rep, optimal_algebra, costs


def normal_form_skew_sym_matrix(mat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the frequencies of the invariant planes of a skew-symmetric matrix. We do so by computing the Schur
    decomposition of the matrix, which gives us approximately a 2x2 block-diagonal matrix. These blocks are skew-sym,
    hence associated to a value, which is the frequency of the invariant plane. We then permute the entries of the
    matrix to make the frequencies increasing, and to have the positive frequencies at the top right of each block.

    Args:
        mat (np.ndarray): A real skew-symmetric matrix of shape (n, n).

    Returns:
        weights (np.ndarray): The sorted, normalized frequencies (length n/2).
        change_of_basis (np.ndarray): The orthogonal matrix that block-diagonalizes the input matrix.
    """
    # Projects on the skew-symmetric matrices
    mat = (mat - mat.T) / 2

    # Computes the Schur decomposition as
    # change_of_basis @ block_diag @ change_of_basis.block_diag = mat.
    block_diag, change_of_basis = scipy.linalg.schur(mat)

    # Extracts the weights of the invariant planes
    # (top right entry of each block)
    ambient_dim = np.shape(mat)[0]
    weights = [block_diag[i, i + 1] for i in range(0, ambient_dim, 2)]
    weights /= np.linalg.norm(weights)

    # Finds location of negative weights, and permute with the positive ones
    index = np.where(weights < 0)[0]
    for i in index:
        weights[i] = -weights[i]
        change_of_basis[:, 2 * i + 1] *= -1

    # Sorts frequencies (permute the matrix)
    index = np.argsort(weights)
    weights = weights[index]
    change_of_basis = change_of_basis[
        :, np.ravel(np.column_stack((2 * index, 2 * index + 1)))
    ]

    return weights, change_of_basis


def normal_form_skew_symmetric_matrices(
    matrices: list[np.ndarray], *, rng=None, initial_point=None,
    optimizer_options=None, runs=None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Computes the weights of a joint normal form of a list of skew-symmetric matrices. This is done by optimizing over
    the changes of basis (matrices in SO(n)) to cancel the entries that are not the antidiagonal terms of the 2x2
    diagonal blocks of the matrices. We assume that the matrices commute, otherwise such a joint normal form does not
    exist.

    Args:
        matrices (list[np.ndarray]): List of real skew-symmetric matrices of shape (n, n).

    Returns:
        weights (np.ndarray): Tuple of tuples containing the sorted frequencies for each matrix.
        change_of_basis (np.ndarray): The orthogonal matrix that block-diagonalizes the input matrices.
    """
    ambient_dim = np.shape(matrices[0])[0]
    m = int(np.shape(matrices[0])[0] / 2)
    # Defines cost function
    manifold = pymanopt.manifolds.SpecialOrthogonalGroup(ambient_dim, k=1)

    @pymanopt.function.autograd(manifold)
    def _cost_function(orth):
        # Conjugate canonical matrices
        matrices_conjugate = [orth @ A @ orth.T for A in matrices]
        # Get diagonal-block terms.
        entries = [
            matrices_conjugate[i][2 * j + k, 2 * j + 1 - k]
            for i in range(len(matrices))
            for j in range(m)
            for k in (0, 1)
        ]
        entries = np.array(entries)

        # Computes norm of not block-diagonal terms.
        differences = np.array([np.trace(A @ A.T) for A in matrices_conjugate])
        norm = np.sum(differences @ differences.T) - np.sum(entries @ entries.T)

        return norm

    # Runs optimization
    problem = pymanopt.Problem(manifold, _cost_function)
    result = _run_optimizer(problem, ambient_dim, rng, initial_point, optimizer_options)
    if runs is not None:
        runs["joint_normal_form"] = _run_summary(result, optimizer_options)
    change_of_basis = result.point
    # Transforms in normal form
    matrices_normal_form = [change_of_basis @ A @ change_of_basis.T for A in matrices]
    weights = tuple(
        [
            tuple([matrices_normal_form[i][2 * j + 0, 2 * j + 1] for j in range(m)])
            for i in range(len(matrices))
        ]
    )

    return weights, change_of_basis.T
