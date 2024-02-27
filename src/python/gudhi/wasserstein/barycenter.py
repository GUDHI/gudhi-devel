# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Theo Lacombe
#
# Copyright (C) 2019 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification


import numpy as np
import warnings

from gudhi.wasserstein import wasserstein_distance


def _mean(x, m):
    '''
    :param x: a list of 2D-points, off diagonal, x_0... x_{k-1}
    :param m: total amount of points taken into account, that is we have (m-k) copies of diagonal
    :returns: the weighted mean of x with (m-k) copies of the diagonal
    '''
    k = len(x)
    if k > 0:
        w = np.mean(x, axis=0)
        w_delta = (w[0] + w[1]) / 2 * np.ones(2)
        return (k * w + (m-k) * w_delta) / m
    else:
        return np.array([0, 0])


def lagrangian_barycenter(pdiagset, init=None, verbose=False):
    '''
    :param pdiagset: a list of ``numpy.array`` of shape `(n x 2)` (`n` can variate), encoding a set of persistence
        diagrams with only finite coordinates.
    :param init: The initial value for barycenter estimate.
        If ``None``, init is made on a random diagram from the dataset.
        Otherwise, it can be an ``int`` (then initialization is made on ``pdiagset[init]``)
        or a `(n x 2)` ``numpy.array`` encoding a persistence diagram with `n` points.
    :type init: ``int``, or (n x 2) ``np.array``
    :param verbose: if ``True``, returns additional information about the barycenter.
    :type verbose: boolean
    :returns: If not verbose (default), a ``numpy.array`` encoding the barycenter estimate of pdiagset
        (local minimum of the energy function).
        If ``pdiagset`` is empty, returns ``None``.
        If verbose, returns a couple ``(Y, log)`` where ``Y`` is the barycenter estimate,
        and ``log`` is a ``dict`` that contains additional information:

        - `"groupings"`, a list of list of pairs ``(i,j)``. Namely, ``G[k] = [...(i, j)...]``, where ``(i,j)``
          indicates that `pdiagset[k][i]` is matched to ``Y[j]`` if ``i = -1`` or ``j = -1``, it means they represent
          the diagonal.

        - `"energy"`, ``float`` representing the Frechet energy value obtained. It is the mean of squared distances of
          observations to the output.

        - `"nb_iter"`, ``int`` number of iterations performed before convergence of the algorithm.
    '''
    X = pdiagset  # to shorten notations, not a copy
    m = len(X)  # number of diagrams we are averaging
    if m == 0:
        warnings.warn("Computing barycenter of empty diag set. Returns None.")
        return None

    # Initialisation of barycenter
    if init is None:
        i0 = np.random.randint(m)  # Index of first state for the barycenter
        Y = X[i0].copy()
    else:
        if isinstance(init, int):
            Y = X[init].copy()
        else:
            Y = init.copy()

    nb_iter = 0

    converged = False  # stopping criterion
    while not converged:
        nb_iter += 1
        K = len(Y)  # current nb of points in Y (some might be on diagonal)
        G = np.full((K, m), -1, dtype=int)  # will store for each j, the (index)
                              # point matched in each other diagram
                              #(might be the diagonal).
                              # that is G[j, i] = k <=> y_j is matched to
                              # x_k in the diagram i-th diagram X[i]
        updated_points = np.zeros((K, 2))  # will store the new positions of
                                           # the points of Y.
                                           # If points disappear, there thrown
                                           # on [0,0] by default.
        new_created_points = []  # will store potential new points.

        # Step 1 : compute optimal matching (Y, X_i) for each X_i
        #          and create new points in Y if needed
        for i in range(m):
            _, indices = wasserstein_distance(Y, X[i], matching=True, order=2., internal_p=2.)
            for y_j, x_i_j in indices:
                if y_j >= 0:  # we matched an off diagonal point to x_i_j...
                    if x_i_j >= 0:  # ...which is also an off-diagonal point.
                        G[y_j, i] = x_i_j
                    else:  # ...which is a diagonal point
                        G[y_j, i] = -1  # -1 stands for the diagonal (mask)
                else:  # We matched a diagonal point to x_i_j...
                    if x_i_j >= 0:  # which is a off-diag point !
                                                # need to create new point in Y
                        new_y = _mean(np.array([X[i][x_i_j]]), m)
                        # Average this point with (m-1) copies of Delta
                        new_created_points.append(new_y)

        # Step 2 : Update current point position thanks to groupings computed
        to_delete = []
        for j in range(K):
            matched_points = [X[i][G[j, i]] for i in range(m) if G[j, i] > -1]
            new_y_j = _mean(matched_points, m)
            if not np.array_equal(new_y_j, np.array([0,0])):
                updated_points[j] = new_y_j
            else: # this points is no longer of any use.
                to_delete.append(j)
        # we remove the point to be deleted now.
        updated_points = np.delete(updated_points, to_delete, axis=0)

        # we cannot converge if there have been new created points.
        if new_created_points:
            Y = np.concatenate((updated_points, new_created_points))
        else:
            # Step 3 : we check convergence
            if np.array_equal(updated_points, Y):
                converged = True
            Y = updated_points


    if verbose:
        groupings = []
        energy = 0
        log = {}
        for i in range(m):
            cost, edges = wasserstein_distance(Y, X[i], matching=True, order=2., internal_p=2.)
            groupings.append(edges)
            energy += cost
            log["groupings"] = groupings
        energy = energy/m
        log["energy"] = energy
        log["nb_iter"] = nb_iter

        return Y, log
    else:
        return Y
