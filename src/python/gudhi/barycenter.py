import ot
import numpy as np
import scipy.spatial.distance as sc


# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Theo Lacombe
#
# Copyright (C) 2019 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification


def _proj_on_diag(w):
    '''
        Util function to project a point on the diag.
    '''
    return np.array([(w[0] + w[1])/2 , (w[0] + w[1])/2])


def _proj_on_diag_array(X):
    '''
    :param X: (n x 2) array encoding the points of a persistent diagram.
    :returns: (n x 2) array encoding the (respective orthogonal) projections of the points onto the diagonal
    '''
    Z = (X[:,0] + X[:,1]) / 2.
    return np.array([Z , Z]).T


def _build_dist_matrix(X, Y, p=2., q=2.):
    '''
    :param X: (n x 2) numpy.array encoding the (points of the) first diagram.
    :param Y: (m x 2) numpy.array encoding the second diagram.
    :param q: Ground metric (i.e. norm l_q).
    :param p: exponent for the Wasserstein metric.
    :returns: (n+1) x (m+1) np.array encoding the cost matrix C. 
                For 1 <= i <= n, 1 <= j <= m, C[i,j] encodes the distance between X[i] and Y[j], while C[i, m+1] (resp. C[n+1, j]) encodes the distance (to the p) between X[i] (resp Y[j]) and its orthogonal proj onto the diagonal.
                note also that C[n+1, m+1] = 0  (it costs nothing to move from the diagonal to the diagonal).
                Note that for lagrangian_barycenter, one must use p=q=2. 
    '''
    Xdiag = _proj_on_diag_array(X)
    Ydiag = _proj_on_diag_array(Y)
    if np.isinf(q):
        C = sc.cdist(X, Y, metric='chebyshev')**p
        Cxd = np.linalg.norm(X - Xdiag, ord=q, axis=1)**p
        Cdy = np.linalg.norm(Y - Ydiag, ord=q, axis=1)**p
    else:
        C = sc.cdist(X,Y, metric='minkowski', p=q)**p
        Cxd = np.linalg.norm(X - Xdiag, ord=q, axis=1)**p
        Cdy = np.linalg.norm(Y - Ydiag, ord=q, axis=1)**p
    Cf = np.hstack((C, Cxd[:,None]))
    Cdy = np.append(Cdy, 0)

    Cf = np.vstack((Cf, Cdy[None,:]))

    return Cf


def _optimal_matching(X, Y):
    """
    :param X: numpy.array of size (n x 2)
    :param Y: numpy.array of size (m x 2)
    :returns: numpy.array of shape (k x 2) encoding the list of edges in the optimal matching. 
                That is, [[(i, j) ...], where (i,j) indicates that X[i] is matched to Y[j]
                if i > len(X) or j > len(Y), it means they represent the diagonal.
              
    """

    n = len(X)
    m = len(Y)
    if X.size == 0: # X is empty
        if Y.size == 0: # Y is empty
            return np.array([[0,0]])  # the diagonal is matched to the diagonal and that's it...
        else:
            return np.column_stack([np.zeros(m+1, dtype=int), np.arange(m+1, dtype=int)]) # TO BE CORRECTED
    elif Y.size == 0: # X is not empty but Y is empty
        return np.column_stack([np.zeros(n+1, dtype=int), np.arange(n+1, dtype=int)])  # TO BE CORRECTED

    # we know X, Y are not empty diags now
    M = _build_dist_matrix(X, Y)

    a = np.full(n+1, 1. / (n + m) )  # weight vector of the input diagram. Uniform here.
    a[-1] = a[-1] * m                # normalized so that we have a probability measure, required by POT
    b = np.full(m+1, 1. / (n + m) )  # weight vector of the input diagram. Uniform here.
    b[-1] = b[-1] * n                # so that we have a probability measure, required by POT
    P = ot.emd(a=a, b=b, M=M)*(n+m)
    # Note : it seems POT return a permutation matrix in this situation,
    # ...guarantee...?
    # It should be enough to check that the algorithm only iterates on vertices of the transportation polytope.
    P[P < 0.5] = 0  # dirty trick to avoid some numerical issues... to be improved.
    # return the list of (i,j) such that P[i,j] > 0, i.e. x_i is matched to y_j (should it be the diag). 
    res = np.nonzero(P)
    return np.column_stack(res)


def _mean(x, m):
    """
    :param x: a list of 2D-points, off diagonal, x_0... x_{k-1}
    :param m: total amount of points taken into account, that is we have (m-k) copies of diagonal
    :returns: the weighted mean of x with (m-k) copies of the diagonal
    """
    k = len(x)
    if k > 0:
        w = np.mean(x, axis=0)
        w_delta = _proj_on_diag(w)
        return (k * w + (m-k) * w_delta) / m
    else:
        return np.array([0, 0])


def lagrangian_barycenter(pdiagset, init=None, verbose=False):
    """
        Compute the estimated barycenter computed with the algorithm provided
        by Turner et al (2014).
        It is a local minima of the corresponding Frechet function.
    :param pdiagset: a list of size N containing numpy.array of shape (n x 2) 
                        (n can variate), encoding a set of 
                        persistence diagrams with only finite coordinates. 
    :param init: The initial value for barycenter estimate. 
                    If None, init is made on a random diagram from the dataset. 
                    Otherwise, it must be an int (then we init with diagset[init])
                    or a (n x 2) numpy.array enconding a persistence diagram with n points.
    :param verbose: if True, returns additional information about the
                        barycenters (assignment and energy).
    :returns: If not verbose (default), a numpy.array encoding
                    the barycenter estimate (local minima of the energy function). 
                    If verbose, returns a triplet (Y, a, e)
                    where Y is the barycenter estimate, a is the assignments between the
                    points of Y and thoses of the diagrams, 
                    and e is the energy value reached by the estimate.
    """
    X = pdiagset  # to shorten notations, not a copy
    m = len(X)  # number of diagrams we are averaging
    if m == 0:
        print("Warning: computing barycenter of empty diag set. Returns None")
        return None

    nb_off_diag = np.array([len(X_i) for X_i in X])  # store the number of off-diagonal point for each of the X_i

    # Initialisation of barycenter
    if init is None:
        i0 = np.random.randint(m)  # Index of first state for the barycenter
        Y = X[i0].copy() #copy() ensure that we do not modify X[i0]
    else:
        if type(init)==int:
            Y = X[init].copy()
        else:
            Y = init.copy()

    converged = False  # stoping criterion
    while not converged:
        K = len(Y)  # current nb of points in Y (some might be on diagonal)
        G = np.zeros((K, m), dtype=int)-1  # will store for each j, the (index) point matched in each other diagram (might be the diagonal). 
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
                    if x_i_j < nb_off_diag[i]:  # ...which is also an off-diagonal point
                        G[y_j, i] = x_i_j
                    else:  # ...which is a diagonal point
                        G[y_j, i] = -1  # -1 stands for the diagonal (mask)
                else:  # We matched a diagonal point to x_i_j...
                    if x_i_j < nb_off_diag[i]:  # which is a off-diag point ! so we need to create a new point in Y
                        new_y = _mean(np.array([X[i][x_i_j]]), m)  # Average this point with (m-1) copies of Delta
                        new_created_points.append(new_y)

        # Step 2 : Update current point position thanks to the groupings computed

        to_delete = []
        for j in range(K):
            matched_points = [X[i][G[j, i]] for i in range(m) if G[j, i] > -1]
            new_y_j = _mean(matched_points, m)
            if not np.array_equal(new_y_j, np.array([0,0])):
                updated_points[j] = new_y_j 
            else: # this points is no longer of any use.
                to_delete.append(j)
        # we remove the point to be deleted now.
        updated_points = np.delete(updated_points, to_delete, axis=0)  # cannot be done in-place.


        if new_created_points: # we cannot converge if there have been new created points.
            Y = np.concatenate((updated_points, new_created_points))
        else:
            # Step 3 : we check convergence
            if np.array_equal(updated_points, Y):
                converged = True 
            Y = updated_points


    if verbose:
        matchings = []
        #energy = 0
        n_y = len(Y)
        for i in range(m):
            edges = _optimal_matching(Y, X[i])
            matchings.append([x_i_j for (y_j, x_i_j) in enumerate(edges) if y_j < n_y])
            # energy += sum([M[i,j] for i,j in enumerate(edges)])

        # energy = energy/m
        return Y, matchings   #, energy
    else:
        return Y

def _plot_barycenter(X, Y, matchings):
    """
    :param X: list of persistence diagrams.
    :param Y: numpy.array of (n x 2). Aims to be an estimate of the barycenter
        returned by lagrangian_barycenter(X, verbose=True).
    :param matchings: list of lists, such that L[k][i] = j if and only if
        the i-th point of the barycenter is grouped with the j-th point of the k-th
        diagram.
    """
    # import matplotlib now to avoid useless dependancies 
    
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # n_y = len(Y.points)
    for i in range(len(X)):
        indices = matchings[i]
        n_i = len(X[i])

        for (y_j, x_i_j) in indices:
            y = Y[y_j]
            if y[0] != y[1]:
                if x_i_j < n_i:  # not mapped with the diag
                    x = X[i][x_i_j]
                else:  # y_j is matched to the diagonal
                    x = _proj_on_diag(y)
                ax.plot([y[0], x[0]], [y[1], x[1]], c='black',
                        linestyle="dashed")
    
    ax.scatter(Y[:,0], Y[:,1], color='purple', marker='d', zorder=2)

    for X_i in X:
        if X_i.size > 0:
            ax.scatter(X_i[:,0], X_i[:,1], marker ='o', zorder=2)

    shift = 0.1  # for improved rendering
    try:
        xmin = np.min(np.array([np.min(x[:,0]) for x in X if len(x) > 0]) - shift)
        xmax = np.max(np.array([np.max(x[:,0]) for x in X if len(x) > 0]) + shift)
        ymin = np.min(np.array([np.max(x[:,1]) for x in X if len(x) > 0]) - shift)
        ymax = np.max(np.array([np.max(x[:,1]) for x in X if len(x) > 0]) + shift)
    except ValueError: # to handle the pecular case where we only average empty diagrams.
        xmin, xmax, ymin, ymax = 0, 1, 0, 1
    themin = min(xmin, ymin)
    themax = max(xmax, ymax)
    ax.set_xlim(themin, themax)
    ax.set_ylim(themin, themax)
    ax.add_patch(Polygon([[themin,themin], [themax,themin], [themax,themax]], fill=True, color='lightgrey'))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal', adjustable='box')
    ax.set_title("Estimated barycenter")

    plt.show()


def _test_perf():
    nb_repeat = 10
    nb_points_in_dgm = [5, 10, 20, 50, 100]
    nb_dmg = [3, 5, 10, 20]

    from time import time
    for m in nb_dmg:
        for n in nb_points_in_dgm:
            tstart = time()
            for _ in range(nb_repeat):
                X = [np.random.rand(n, 2) for _ in range(m)]
                for diag in X:
                    #enforce having diagrams
                    diag[:,1] = diag[:,1] + diag[:,0]
                _ = lagrangian_barycenter(X)
            tend = time()
            print("Computation of barycenter in %s sec, with k = %s diags and n = %s points per diag."%(np.round((tend - tstart)/nb_repeat, 2), m, n))
            print("********************")


def _sanity_check(verbose):
    #dg1 = np.array([[0.2, 0.5]])
    #dg2 = np.array([[0.2, 0.7], [0.73, 0.88]])
    #dg3 = np.array([[0.3, 0.6], [0.7, 0.85], [0.2, 0.3]])
    #X = [dg1, dg2, dg3]
    #Y, a = lagrangian_barycenter(X, verbose=verbose)
    #_plot_barycenter(X, Y, a)

    #dg1 = np.array([[0.2, 0.5]])
    #dg2 = np.array([]) # The empty diagram
    #dg3 = np.array([[0.4, 0.8]])
    #X = [dg1, dg2, dg3]
    #Y, a = lagrangian_barycenter(X, verbose=verbose)
    #_plot_barycenter(X, Y, a)
 
    #dg1 = np.array([])
    #dg2 = np.array([]) # The empty diagram
    #dg3 = np.array([])
    #X = [dg1, dg2, dg3]
    #Y, a = lagrangian_barycenter(X, verbose=verbose)
    #_plot_barycenter(X, Y, a)
    #print(Y)
 
    dg1 = np.array([[0.1, 0.12], [0.21, 0.7], [0.4, 0.5], [0.3, 0.4], [0.35, 0.7], [0.5, 0.55], [0.32, 0.42], [0.1, 0.4], [0.2, 0.4]])
    dg2 = np.array([[0.09, 0.11], [0.3, 0.43], [0.5, 0.61], [0.3, 0.7], [0.42, 0.5], [0.35, 0.41], [0.74, 0.9], [0.5, 0.95], [0.35, 0.45], [0.13, 0.48], [0.32, 0.45]])
    dg3 = np.array([[0.1, 0.15], [0.1, 0.7], [0.2, 0.22], [0.55, 0.84], [0.11, 0.91], [0.61, 0.75], [0.33, 0.46], [0.12, 0.41], [0.32, 0.48]])
    X = [dg3]
    Y, a = lagrangian_barycenter(X, verbose=verbose)
    _plot_barycenter(X, Y, a)
    print(Y)
 
    
    #dg1 = np.array([[0.2, 0.5]])
    #dg2 = np.array([[0.2, 0.7]])
    #dg3 = np.array([[0.3, 0.6], [0.7, 0.8], [0.2, 0.3]])
    #dg4 = np.array([])
    #
    #bary, a = lagrangian_barycenter(pdiagset=[dg1, dg2, dg3, dg4],init=3, verbose=True)
    #_plot_barycenter([dg1, dg2, dg3, dg4], bary, a)
    #message = "Wasserstein barycenter estimated:"    
    #print(message)
    #print(bary)

if __name__=="__main__":
    _sanity_check(verbose = True)
    #_test_perf()
