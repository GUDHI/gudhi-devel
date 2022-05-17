# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Theo Lacombe
#
# Copyright (C) 2021 Université Gustave Eiffel
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification


import numpy as np
from _utils import _dist_to_diag, _build_dist_matrix
import warnings


def _get_cells(X, c, withdiag, internal_p):
    """
    Given a point cloud `X` and a codebook (set of centroids) `c`,
    assign to each point in X its nearest centroid in c.
    If withdiag is True, points which are nearest to the diagonal than any centroid are assigned to an
    additional cell.

    Parameters
    ----------
    :param X: Current set of points, of size (n x d) (n : number of points, and d=2 for persistence diagrams).
    :type X: `numpy.array`
    :param c: Current set of centroids (codebook), of size (k x d). (k centroids in dimension d, d=2 for PDs).
                It is assumed that k >= 1.
    :param withdiag: Indicate if we consider the diagonal into account (True when working with PDs).
    :type withdiag: bool
    :param internal_p: The ground metric used to compute distances between points (norm_p).

    :returns: list `cells` of affectation, such that cells[j] = points in `X` matched to `c[j]`,
                with convention that `cells[k]` represents the diagonal whenever `withdiags=True`.

        .. note: if the input X is empty, returns the empty list. Should not happen as _get_batches remove
                    batches of empty points.

    """
    k = c.shape[0]

    if X.shape[0] == 0:
        warnings.warn("Input point cloud is empty, Cells affectation is a list of empty lists.")
        return [[]] * (k+withdiag)  # k or (k+1) depending on additional centroid for the diagonal.

    if X.shape[1] != 2:
        raise ValueError("Input batch must be of shape (n x 2), not (%s,%s)" %(X.shape))

    M = _build_dist_matrix(X, c, order=2, internal_p=internal_p)  # Note: Order is useless here
                                                                  # we just need nearest neighbors.

    if withdiag:
        a = np.argmin(M[:-1, :], axis=1)
    else:
        a = np.argmin(M[:-1, :-1], axis=1)

    cells = [X[a == j] for j in range(k)]

    if withdiag:
        cells.append(X[a == k])  # this is the (k+1)-th centroid corresponding to the diagonal

    return cells


def _from_batch(pdiagset, batches_indices):
    """
    Given a pdiagset and a list of indices, compute the "empirical expected persistence diagram"
    (in practice, concatenates) on the corresponding batch of diagrams.
    Note that we may encounter empty diagrams, hence the specific check .ndim==2.
        If the resulting list is empty (all diagrams are empty), returns empty array.
    """
    list_of_non_empty_diags = [pdiagset[i] for i in batches_indices if pdiagset[i].ndim==2]
    if list_of_non_empty_diags: # There are some non-empty diagrams
        X_batch = np.concatenate(list_of_non_empty_diags)
        return X_batch
    else:
        return np.array([])


def _init_c(pdiagset, k, internal_p=2):
    """
    A naive heuristic to initialize a codebook: we take the k points with largest distances to the diagonal
        in the first diagram of the list.

    Parameters
    ----------
    :param list_diags: The list of diagrams we want to quantize. Assume it is not empty (checked by quantization).
    :param k: number of centroids in the codebook.
    :param internal_p: Internal ground metric. Should not play a role specific role.

    :returns: an initial codebook c0 with `k` entries.
    -------

    """
    dgm = pdiagset[0]
    w = _dist_to_diag(dgm, internal_p)
    s = np.argsort(w)
    c0 = dgm[s[-k:]]  # Elements in dgm with larger persistence are last
    return c0


#####################
###  Main method  ###
#####################

def quantization(pdiagset, k=2, init=None, batch_size=1, order=2., internal_p=2.):
    """
    This quantization algorithm takes a list of diagrams ``pdiagset``, an integer ``k``
    (or an initial codebook guess ``c0``)
    and returns a codebook of ``k`` centroids that aims at minimizing the distance between the codebook and the
    expected persistence diagram underneath the diagrams in ``pdiagset``.
    The codebook is iteratively updated by going through (batches of) persistence diagrams.

    :param pdiagset: a set of persistence diagrams that we want to quantize.
    :type pdiagset: list of n x 2 numpy arrays
    :param k: number of centroids. Default is ``2``. A naive heuristic is used to initialize the codebook.
              Not used if an initial codebook ``c0`` has been provided.
    :type k: ``int``
    :param init: If provided, overwrite the provided ``k`` (which is now ``c0.shape[0]``), and is used as an
                initialization for the algorithm. Default is ``None``.
    :type init: (n x 2) ``numpy.array``, or ``None``.
    :param batch_size: Size of batches used during the online exploration of the ``pdiagset``.
                        Default is ``1``.
    :type batch_size: ``int``
    :param order: Order of the Wasserstein distance we want to minimize.
                  Only 2 is implemented for now. Default is ``2.``.
    :type order: ``float``
    :param internal_p: Ground metric to assess centroid affectation. Default is ``2.``.
    :type internal_p: ``float``

    :returns: The final codebook obtained after going through the all pdiagset.

    .. note:: The exact algorithm presented in the reference paper
              requires more technical considerations in order to prove theoretical convergence rates.
              A theoretically motivated value (asymptotically optimal) for the ``batch_size`` is
              ``int(log(len(pdiagset)))``.
              However, numerical experiments suggest that this simplified version yields substantially similar results.
    """

    # We take the diagonal into account in the following.
    # This could be turned into a parameter of the algorithm (so that it would encompass the standard
    # "online-Lloyd algorithm"), but I thought this may be confusing.
    withdiag = True

    if not pdiagset:
        raise ValueError("Input pdiagset is empty.")
    # If no codebook has been provided, we propose an initialization.
    if init is None:
        init = _init_c(pdiagset, k)
    # Number of centroids in the codebook (updated in the case c0 was provided by the user).
    k = init.shape[0]

    # Variable that will store the centroid evolution.
    c_current = init.copy()

    # number of diagrams
    n = len(pdiagset)

    # Building the batches (size batch_size, except possibly for the last one): a list of the form
    # [[0,1,..., batch_size-1],[batch_size,..., 2batch_size-1],...,[remaining diags]]
    # Note : we have to add (n % batch_size != 0) to increase the number of batches by 1 if needed.
    nb_step = n // batch_size + (n % batch_size != 0)
    batches = np.array_split(np.arange(0, n, dtype=int), nb_step)

    # We now iterate through the diagrams
    for t in range(nb_step):
        # Get the empirical expected persistence diagram on the current batch
        X_bar_t = _from_batch(pdiagset, batches[t])
        # We compute the corresponding cells
        cells = _get_cells(X_bar_t, c_current, withdiag=withdiag, internal_p=internal_p)

        if order == 2.:
            # We update the position of each centroid.
            for j in range(k):
                lc1 = len(cells[j])
                # The batch may be empty (if c_current[j] does not catch any point)
                # In which case we do not update it.
                if lc1 > 0:
                    # We compute the "gradient", i.e. the optimal descent direction for c_current[j]
                    # When order == 2, it just consists of pushing c_current[j] towards the barycenter
                    # of the cells[j] it defines.
                    grad = np.mean(c_current[j] - cells[j], axis=0)
                    # We then apply the gradient, with a 1/(1+t) factor to guarantee the convergence of this
                    # stochastic-gradient-descent like approach (decreasing learning rate).
                    c_current[j] = c_current[j] - grad / (t + 1)
        else:
            raise NotImplemented('Order = %s is not available yet. Only order=2. is valid' %order)

    return c_current
