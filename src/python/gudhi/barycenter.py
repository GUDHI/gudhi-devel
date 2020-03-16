# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Theo Lacombe
#
# Copyright (C) 2019 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification


import ot
import numpy as np
import scipy.spatial.distance as sc

from gudhi.wasserstein import _build_dist_matrix, _perstot



def _mean(x, m):
    '''
    :param x: a list of 2D-points, off diagonal, x_0... x_{k-1}
    :param m: total amount of points taken into account, 
                that is we have (m-k) copies of diagonal
    :returns: the weighted mean of x with (m-k) copies of the diagonal
    '''
    k = len(x)
    if k > 0:
        w = np.mean(x, axis=0)
        w_delta = (w[0] + w[1]) / 2 * np.ones(2)
        return (k * w + (m-k) * w_delta) / m
    else:
        return np.array([0, 0])


def _optimal_matching(X, Y, withcost=False):
    '''
    :param X: numpy.array of size (n x 2)
    :param Y: numpy.array of size (m x 2)
    :param withcost: returns also the cost corresponding to the optimal matching
    :returns: numpy.array of shape (k x 2) encoding the list of edges 
                                            in the optimal matching. 
                That is, [[i, j] ...], where (i,j) indicates 
                that X[i] is matched to Y[j]
                if i >= len(X) or j >= len(Y), it means they 
                represent the diagonal.
                They will be encoded by -1 afterwards.

    NOTE : this code will be removed for final merge, 
    and wasserstein.optimal_matching will be used instead.
    '''

    n = len(X)
    m = len(Y)
    # Start by handling empty diagrams. Could it be shorten?
    if X.size == 0: # X is empty
        if Y.size == 0: # Y is empty
            res = np.array([[0,0]])  # the diagonal is matched to the diagonal
            if withcost:
                return res, 0
            else:
                return res
        else: # X is empty but not Y
            res = np.array([[0, i] for i in range(m)]) 
            cost = _perstot(Y, order=2, internal_p=2)**2
            if withcost:
                return res, cost
            else:
                return res
    elif Y.size == 0: # X is not empty but Y is empty
        res = np.array([[i,0] for i in range(n)])        
        cost = _perstot(X, order=2, internal_p=2)**2
        if withcost:
            return res, cost
        else:
            return res

    # we know X, Y are not empty diags now
    M = _build_dist_matrix(X, Y, order=2, internal_p=2)

    a = np.ones(n+1)
    a[-1] = m
    b = np.ones(m+1)
    b[-1] = n
    P = ot.emd(a=a, b=b, M=M) 
    # Note : it seems POT returns a permutation matrix in this situation, 
    # ie a vertex of the constraint set (generically true).
    if withcost:
        cost = np.sum(np.multiply(P, M))
    P[P < 0.5] = 0  # dirty trick to avoid some numerical issues... to improve.
    res = np.argwhere(P)

    # return the list of (i,j) such that P[i,j] > 0, 
    #i.e. x_i is matched to y_j (should it be the diag). 
    if withcost:
        return res, cost
    return res


def lagrangian_barycenter(pdiagset, init=None, verbose=False):
    '''
    :param pdiagset: a list of size m containing numpy.array of shape (n x 2) 
                    (n can variate), encoding a set of 
                    persistence diagrams with only finite coordinates. 
    :param init: The initial value for barycenter estimate. 
                    If None, init is made on a random diagram from the dataset. 
                    Otherwise, it must be an int 
                    (then we init with diagset[init])
                    or a (n x 2) numpy.array enconding 
                    a persistence diagram with n points.
    :param verbose: if True, returns additional information about the
                    barycenter.
    :returns: If not verbose (default), a numpy.array encoding
                    the barycenter estimate of pdiagset
                    (local minima of the energy function). 
                    If pdiagset is empty, returns None.
                    If verbose, returns a couple (Y, log)
                    where Y is the barycenter estimate,
                    and log is a dict that contains additional informations:
                    - groupings, a list of list of pairs (i,j),
                    That is, G[k] = [(i, j) ...], where (i,j) indicates 
                    that X[i] is matched to Y[j]
                    if i = -1 or j = -1, it means they 
                    represent the diagonal.
                    - energy, a float representing the Frechet 
                    energy value obtained,
                    that is the mean of squared distances 
                    of observations to the output.
                    - nb_iter, integer representing the number of iterations 
                    performed before convergence of the algorithm.
    '''
    X = pdiagset  # to shorten notations, not a copy
    m = len(X)  # number of diagrams we are averaging
    if m == 0:
        print("Warning: computing barycenter of empty diag set. Returns None")
        return None
    
    # store the number of off-diagonal point for each of the X_i
    nb_off_diag = np.array([len(X_i) for X_i in X])  
    # Initialisation of barycenter
    if init is None:
        i0 = np.random.randint(m)  # Index of first state for the barycenter
        Y = X[i0].copy() 
    else:
        if type(init)==int:
            Y = X[init].copy()
        else:
            Y = init.copy()

    nb_iter = 0

    converged = False  # stoping criterion
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
            indices = _optimal_matching(Y, X[i])
            for y_j, x_i_j in indices:
                if y_j < K:  # we matched an off diagonal point to x_i_j...
                # ...which is also an off-diagonal point.
                    if x_i_j < nb_off_diag[i]:  
                        G[y_j, i] = x_i_j
                    else:  # ...which is a diagonal point
                        G[y_j, i] = -1  # -1 stands for the diagonal (mask)
                else:  # We matched a diagonal point to x_i_j...
                    if x_i_j < nb_off_diag[i]:  # which is a off-diag point ! 
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
        n_y = len(Y)
        for i in range(m):
            edges, cost = _optimal_matching(Y, X[i], withcost=True)
            n_x = len(X[i])
            G = edges[np.where(edges[:,0]<n_y)]
            idx = np.where(G[:,1] >= n_x)
            G[idx,1] = -1  # -1 will encode the diagonal
            groupings.append(G)
            energy += cost
            log["groupings"] = groupings
        energy = energy/m
        log["energy"] = energy
        log["nb_iter"] = nb_iter

        return Y, log
    else:
        return Y

