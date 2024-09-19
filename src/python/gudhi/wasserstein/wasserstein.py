# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Theo Lacombe
#
# Copyright (C) 2019 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np
import scipy.spatial.distance as sc
import warnings

import ot

# Currently unused, but Théo says it is likely to be used again.
def _proj_on_diag(X):
    '''
    :param X: (n x 2) array encoding the points of a persistent diagram.
    :returns: (n x 2) array encoding the (respective orthogonal) projections of the points onto the diagonal
    '''
    Z = (X[:,0] + X[:,1]) / 2.
    return np.array([Z , Z]).T


def _dist_to_diag(X, internal_p):
    '''
    :param X: (n x 2) array encoding the points of a persistent diagram.
    :param internal_p: Ground metric (i.e. norm L^p).
    :returns: (n) array encoding the (respective orthogonal) distances of the points to the diagonal

    .. note::
        Assumes that the points are above the diagonal.
    '''
    return (X[:, 1] - X[:, 0]) * 2 ** (1.0 / internal_p - 1)


def _build_dist_matrix(X, Y, order, internal_p):
    '''
    :param X: (n x 2) numpy.array encoding the (points of the) first diagram.
    :param Y: (m x 2) numpy.array encoding the second diagram.
    :param order: exponent for the Wasserstein metric.
    :param internal_p: Ground metric (i.e. norm L^p).
    :returns: (n+1) x (m+1) np.array encoding the cost matrix C.
                For 0 <= i < n, 0 <= j < m, C[i,j] encodes the distance between X[i] and Y[j],
                while C[i, m] (resp. C[n, j]) encodes the distance (to the p) between X[i] (resp Y[j])
                and its orthogonal projection onto the diagonal.
                note also that C[n, m] = 0  (it costs nothing to move from the diagonal to the diagonal).
    '''
    Cxd = _dist_to_diag(X, internal_p)**order
    Cdy = _dist_to_diag(Y, internal_p)**order
    if np.isinf(internal_p):
        C = sc.cdist(X,Y, metric='chebyshev')**order
    else:
        C = sc.cdist(X,Y, metric='minkowski', p=internal_p)**order
    Cf = np.hstack((C, Cxd[:,None]))
    Cdy = np.append(Cdy, 0)

    Cf = np.vstack((Cf, Cdy[None,:]))

    return Cf


def _perstot_autodiff(X, order, internal_p):
    '''
    Version of _perstot that works on eagerpy tensors.
    '''
    return _dist_to_diag(X, internal_p).norms.lp(order)


def _perstot(X, order, internal_p, enable_autodiff):
    '''
    :param X: (n x 2) numpy.array (points of a given diagram).
    :param order: exponent for Wasserstein.
    :param internal_p: Ground metric on the (upper-half) plane (i.e. norm L^p in R^2).
    :param enable_autodiff: If X is torch.tensor, tensorflow.Tensor or jax.numpy.ndarray, make the computation
        transparent to automatic differentiation.
    :type enable_autodiff: bool
    :returns: float, the total persistence of the diagram (that is, its distance to the empty diagram).

    .. note::
        Can be +inf if the diagram has an essential part (points with infinite coordinates).
    '''
    if enable_autodiff:
        import eagerpy as ep

        return _perstot_autodiff(ep.astensor(X), order, internal_p).raw
    else:
        return np.linalg.norm(_dist_to_diag(X, internal_p), ord=order)


def _get_essential_parts(a):
    '''
    :param a: (n x 2) numpy.array (point of a diagram)
    :returns: five lists of indices (between 0 and len(a)) accounting for the five types of points with infinite
    coordinates that can occur in a diagram, namely:
        type0 : (-inf, finite)
        type1 : (finite, +inf)
        type2 : (-inf, +inf)
        type3 : (-inf, -inf)
        type4 : (+inf, +inf)
    .. note::
        For instance, a[_get_essential_parts(a)[0]] returns the points in a of coordinates (-inf, x) for some finite x.
        Note also that points with (+inf, -inf) are not handled (points (x,y) in dgm satisfy by assumption (y >= x)).

        Finally, we consider that points with coordinates (-inf,-inf) and (+inf, +inf) belong to the diagonal.
    '''
    if len(a):
        first_coord_finite = np.isfinite(a[:,0])
        second_coord_finite = np.isfinite(a[:,1])
        first_coord_infinite_positive = (a[:,0] == np.inf)
        second_coord_infinite_positive = (a[:,1] == np.inf)
        first_coord_infinite_negative = (a[:,0] == -np.inf)
        second_coord_infinite_negative = (a[:,1] == -np.inf)

        # coord (-inf, x)
        ess_first_type  = np.where(second_coord_finite & first_coord_infinite_negative)[0]
        # coord (x, +inf)
        ess_second_type = np.where(first_coord_finite & second_coord_infinite_positive)[0]
        # coord (-inf, +inf)
        ess_third_type  = np.where(first_coord_infinite_negative & second_coord_infinite_positive)[0]
        # coord (-inf, -inf)
        ess_fourth_type = np.where(first_coord_infinite_negative & second_coord_infinite_negative)[0]
        # coord (+inf, +inf)
        ess_fifth_type  = np.where(first_coord_infinite_positive  & second_coord_infinite_positive)[0]
        return ess_first_type, ess_second_type, ess_third_type, ess_fourth_type, ess_fifth_type
    else:
        return [], [], [], [], []


def _cost_and_match_essential_parts(X, Y, idX, idY, order, axis):
    '''
    :param X: (n x 2) numpy.array (dgm points)
    :param Y: (n x 2) numpy.array (dgm points)
    :param idX: indices to consider for this one dimensional OT problem (in X)
    :param idY: indices to consider for this one dimensional OT problem (in Y)
    :param order: exponent for Wasserstein distance computation
    :param axis: must be 0 or 1, correspond to the coordinate which is finite.
    :returns: cost (float) and match for points with *one* infinite coordinate.

    .. note::
        Assume idX, idY come when calling _handle_essential_parts, thus have same length.
    '''
    u = X[idX, axis]
    v = Y[idY, axis]

    cost = np.sum(np.abs(np.sort(u) - np.sort(v))**(order))  # OT cost in 1D

    sortidX = idX[np.argsort(u)]
    sortidY = idY[np.argsort(v)]
    # We return [i,j] sorted per value
    match = list(zip(sortidX, sortidY))

    return cost, match


def _handle_essential_parts(X, Y, order):
    '''
    :param X: (n x 2) numpy array, first diagram.
    :param Y: (n x 2) numpy array, second diagram.
    :order: Wasserstein order for cost computation.
    :returns: cost and matching due to essential parts. If cost is +inf, matching will be set to None.
    '''
    ess_parts_X = _get_essential_parts(X)
    ess_parts_Y = _get_essential_parts(Y)

    # Treats the case of infinite cost (cardinalities of essential parts differ).
    for u, v in list(zip(ess_parts_X, ess_parts_Y))[:3]: # ignore types 4 and 5 as they belong to the diagonal
        if len(u) != len(v):
            return np.inf, None

    # Now we know each essential part has the same number of points in both diagrams.
    # Handle type 0 and type 1 essential parts (those with one finite coordinates)
    c1, m1 = _cost_and_match_essential_parts(X, Y, ess_parts_X[0], ess_parts_Y[0], axis=1, order=order)
    c2, m2 = _cost_and_match_essential_parts(X, Y, ess_parts_X[1], ess_parts_Y[1], axis=0, order=order)

    c = c1 + c2
    m = m1 + m2

    # Handle type3 (coordinates (-inf,+inf), so we just align points)
    m += list(zip(ess_parts_X[2], ess_parts_Y[2]))

    # Handle type 4 and 5, considered as belonging to the diagonal so matched to (-1) with cost 0.
    for z in ess_parts_X[3:]:
        m += [(u, -1) for u in z] # points in X are matched to -1
    for z in ess_parts_Y[3:]:
        m += [(-1, v) for v in z] # -1 is match to points in Y

    return c, np.array(m)


def _finite_part(X):
    '''
    :param X: (n x 2) numpy array encoding a persistence diagram.
    :returns: The finite part of a diagram `X` (points with finite coordinates).
    '''
    return X[np.where(np.isfinite(X[:,0]) & np.isfinite(X[:,1]))]


def _warn_infty(matching):
    '''
    Handle essential parts with different cardinalities. Warn the user about cost being infinite and (if
    `matching=True`) about the returned matching being `None`.
    '''
    user_warning = 'Cardinality of essential parts differs.'
    if matching:
        warnings.warn(f'{user_warning} Distance (cost) is +inf, and the returned matching is None.')
        return np.inf, None
    else:
        warnings.warn(f'{user_warning} Distance (cost) is +inf.')
        return np.inf


def wasserstein_distance(X, Y, matching=False, order=1., internal_p=np.inf, enable_autodiff=False,
                         keep_essential_parts=True):
    '''
    Compute the Wasserstein distance between persistence diagram using Python Optimal Transport backend.
    Diagrams can contain points with infinity coordinates (essential parts).
    Points with (-inf,-inf) and (+inf,+inf) coordinates are considered as belonging to the diagonal.
    If the distance between two diagrams is +inf (which happens if the cardinalities of essential
    parts differ) and optimal matching is required, it will be set to ``None``.

    :param X: The first diagram.
    :type X: n x 2 numpy.array
    :param Y:  The second diagram.
    :type Y: m x 2 numpy.array
    :param matching: if ``True``, computes and returns the optimal matching between X and Y, encoded as
        a (n x 2) np.array  [...[i,j]...], meaning the i-th point in X is matched to
        the j-th point in Y, with the convention that (-1) represents the diagonal.
    :param order: Wasserstein exponent q (1 <= q < infinity).
    :type order: float
    :param internal_p: Ground metric on the (upper-half) plane (i.e. norm L^p in R^2).
    :type internal_p: float
    :param enable_autodiff: If X and Y are ``torch.tensor`` or ``tensorflow.Tensor``, make the computation
        transparent to automatic differentiation. This requires the package EagerPy and is currently incompatible
        with ``matching=True`` and with ``keep_essential_parts=True``.

        .. note:: This considers the function defined on the coordinates of the off-diagonal finite points of X and Y
            and lets the various frameworks compute its gradient. It never pulls new points from the diagonal.
    :type enable_autodiff: bool
    :param keep_essential_parts: If ``False``, only considers the finite points in the diagrams.
                                 Otherwise, include essential parts in cost and matching computation.
    :type keep_essential_parts: bool
    :returns: The Wasserstein distance of order q (1 <= q < infinity) between persistence diagrams with
              respect to the internal_p-norm as ground metric.
              If matching is set to True, also returns the optimal matching between X and Y.
              If cost is +inf, any matching is optimal and thus it returns `None` instead.
    '''

    # First step: handle empty diagrams
    n = len(X)
    m = len(Y)

    if n == 0:
        if m == 0:
            if not matching:
                # What if enable_autodiff?
                return 0.
            else:
                return 0., np.array([])
        else:
            cost = _perstot(Y, order, internal_p, enable_autodiff)
            if cost == np.inf:
                return _warn_infty(matching)
            else:
                if not matching:
                    return cost
                else:
                    return cost, np.array([[-1, j] for j in range(m)])
    elif m == 0:
        cost = _perstot(X, order, internal_p, enable_autodiff)
        if cost == np.inf:
            return _warn_infty(matching)
        else:
            if not matching:
                return cost
            else:
                return cost, np.array([[i, -1] for i in range(n)])


    # Check essential part and enable autodiff together
    if enable_autodiff and keep_essential_parts:
        warnings.warn('''enable_autodiff=True and keep_essential_parts=True are incompatible together.
                      keep_essential_parts is set to False: only points with finite coordinates are considered
                      in the following.
                      ''')
        keep_essential_parts = False

    # Second step: handle essential parts if needed.
    if keep_essential_parts:
        essential_cost, essential_matching = _handle_essential_parts(X, Y, order=order)
        if (essential_cost == np.inf):
            return _warn_infty(matching)  # Tells the user that cost is infty and matching (if True) is None.
            # avoid computing transport cost between the finite parts if essential parts
            # cardinalities do not match (saves time)
    else:
        essential_cost = 0
        essential_matching = None

    # Now the standard pipeline for finite parts
    if enable_autodiff:
        import eagerpy as ep

        X_orig = ep.astensor(X)
        Y_orig = ep.astensor(Y)
        X = X_orig.numpy()
        Y = Y_orig.numpy()

    # Extract finite points of the diagrams. 
    X, Y = _finite_part(X), _finite_part(Y)
    n = len(X)
    m = len(Y)

    M = _build_dist_matrix(X, Y, order=order, internal_p=internal_p)
    a = np.ones(n+1) # weight vector of the input diagram. Uniform here.
    a[-1] = m
    b = np.ones(m+1) # weight vector of the input diagram. Uniform here.
    b[-1] = n

    if matching:
        assert not enable_autodiff, "matching and enable_autodiff are currently incompatible"
        P = ot.emd(a=a,b=b,M=M, numItermax=2000000)
        ot_cost = np.sum(np.multiply(P,M))
        P[-1, -1] = 0  # Remove matching corresponding to the diagonal
        match = np.argwhere(P)
        # Now we turn to -1 points encoding the diagonal
        match[:,0][match[:,0] >= n] = -1
        match[:,1][match[:,1] >= m] = -1
        # Finally incorporate the essential part matching
        if essential_matching is not None:
            match = np.concatenate([match, essential_matching]) if essential_matching.size else match
        return (ot_cost + essential_cost) ** (1./order) , match

    if enable_autodiff:
        P = ot.emd(a=a, b=b, M=M, numItermax=2000000)
        pairs_X_Y = np.argwhere(P[:-1, :-1])
        pairs_X_diag = np.nonzero(P[:-1, -1])
        pairs_Y_diag = np.nonzero(P[-1, :-1])
        dists = []
        # empty arrays are not handled properly by the helpers, so we avoid calling them
        if len(pairs_X_Y):
            dists.append((Y_orig[pairs_X_Y[:, 1]] -
                          X_orig[pairs_X_Y[:, 0]]).norms.lp(internal_p, axis=-1).norms.lp(order))
        if len(pairs_X_diag[0]):
            dists.append(_perstot_autodiff(X_orig[pairs_X_diag], order, internal_p))
        if len(pairs_Y_diag[0]):
            dists.append(_perstot_autodiff(Y_orig[pairs_Y_diag], order, internal_p))
        dists = [dist.reshape(1) for dist in dists]
        return ep.concatenate(dists).norms.lp(order).raw
        # We can also concatenate the 3 vectors to compute just one norm.

    # Comptuation of the ot cost using the ot.emd2 library.
    # Note: it is the Wasserstein distance to the power q.
    # The default numItermax=100000 is not sufficient for some examples with 5000 points, what is a good value?
    ot_cost = ot.emd2(a, b, M, numItermax=2000000)

    return (ot_cost + essential_cost) ** (1./order)
