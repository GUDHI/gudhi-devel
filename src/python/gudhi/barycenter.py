import ot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def _proj_on_diag(x):
    return np.array([(x[0] + x[1]) / 2, (x[0] + x[1]) / 2])


def _norm2(x, y):
    return (y[0] - x[0])**2 + (y[1] - x[1])**2


def _norm_inf(x, y):
    return np.max(np.abs(y[0] - x[0]), np.abs(y[1] - x[1]))


def _cost_matrix(X, Y):
    """
    :param X: (n x 2) numpy.array encoding the first diagram
    :param Y: (m x 2) numpy.array encoding the second diagram
    :return: The cost matrix with size (k x k) where k = |d_1| + |d_2| in order to encode matching to diagonal
    """
    n, m = len(X), len(Y)
    k = n + m
    M = np.zeros((k, k))
    for i in range(n):  # go throught X points
        x_i = X[i]
        p_x_i = _proj_on_diag(x_i)  # proj of x_i on the diagonal
        dist_x_delta = _norm2(x_i, p_x_i)  # distance to the diagonal regarding the ground norm
        for j in range(m):  # go throught d_2 points
            y_j = Y[j]
            p_y_j = _proj_on_diag(y_j)
            M[i, j] = _norm2(x_i, y_j)
            dist_y_delta = _norm2(y_j, p_y_j)
            for it in range(m):
                M[n + it, j] = dist_y_delta
        for it in range(n):
            M[i, m + it] = dist_x_delta

    return M


def _optimal_matching(M):
    n = len(M)
    # if input weights are empty lists, pot treat the uniform assignement problem and returns a bistochastic matrix (up to *n).
    P = ot.emd(a=[], b=[], M=M) * n
    # return the list of indices j such that L[i] = j iff P[i,j] = 1
    return np.nonzero(P)[1]


def _mean(x, m):
    """
    :param x: a list of 2D-points, of diagonal, x_0... x_{k-1}
    :param m: total amount of points taken into account, that is we have (m-k) copies of diagonal
    :returns: the weighted mean of x with (m-k) copies of Delta taken into account (defined by mukherjee etc.)
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
        Compute the estimated barycenter computed with the Hungarian algorithm provided by Mukherjee et al
                    It is a local minima of the corresponding Frechet function.
                    It exactly belongs to the persistence diagram space (because all computations are made on it).
        :param pdiagset: a list of size N containing numpy.array of shape (n x
                            2) (n can variate), encoding a set of persistence diagrams with only finite
                            coordinates. 
        :param init: The initial value for barycenter estimate. If None, init is made on a random diagram from the dataset. Otherwise, it must be a (n x 2) numpy.array enconding a persistence diagram with n points.
        :returns: If not verbose (default), the barycenter estimate (local minima of the energy function). If verbose, returns a triplet (Y, a, e) where Y is the barycenter estimate, a is the assignments between the points of Y and thoses of the diagrams, and e is the energy value reached by the estimate.
    """
    m = len(pdiagset)  # number of diagrams we are averaging
    X = pdiagset  # to shorten notations
    nb_off_diag = np.array([len(X_i) for X_i in X])  # store the number of off-diagonal point for each of the X_i

    # Initialisation of barycenter
    if init is None:
        i0 = np.random.randint(m)  # Index of first state for the barycenter
        Y = X[i0].copy()
    else:
        Y = init.copy()

    not_converged = True  # stoping criterion
    while not_converged:
        K = len(Y)  # current nb of points in Y (some might be on diagonal)
        G = np.zeros((K, m))  # will store for each j, the (index) point matched in each other diagram (might be the diagonal). 
        updated_points = np.zeros((K, 2))  # will store the new positions of the points of Y
        new_created_points = []  # will store eventual new points.

        # Step 1 : compute optimal matching (Y, X_i) for each X_i
        for i in range(m):
            M = _cost_matrix(Y, X[i])
            indices = _optimal_matching(M)
            for y_j, x_i_j in enumerate(indices):
                if y_j < K:  # we matched an off diagonal point to x_i_j...
                    if x_i_j < nb_off_diag[i]:  # ...which is also an off-diagonal point
                        G[y_j, i] = x_i_j
                    else:  # ...which is a diagonal point
                        G[y_j, i] = -1  # -1 stands for the diagonal (mask)
                else:  # We matched a diagonal point to x_i_j...
                    if x_i_j < nb_off_diag[i]:  # which is a off-diag point ! so we need to create a new point in Y
                        new_y = _mean(np.array([X[i][x_i_j]]), m)  # Average this point with (m-1) copies of Delta
                        new_created_points.append(new_y)

        # Step 2 : Compute new points (mean)
        for j in range(K):
            matched_points = [X[i][int(G[j, i])] for i in range(m) if G[j, i] > -1]
            updated_points[j] = _mean(matched_points, m)

        if new_created_points:
            Y = np.concatenate((updated_points, new_created_points))
        else:
            Y = updated_points

        # Step 3 : we update our estimation of the barycenter
        if len(new_created_points) == 0 and np.array_equal(updated_points, Y):
            not_converged = False

    if verbose:
        matchings = []
        energy = 0
        n_y = len(Y)
        for i in range(m):
            M = _cost_matrix(Y, X[i])
            edges = _optimal_matching(M)
            matchings.append([x_i_j for (y_j, x_i_j) in enumerate(edges) if y_j < n_y])
            #energy += total_cost

        #energy /= m
        _plot_barycenter(X, Y, matchings)
        plt.show()
        return Y, matchings, energy
    else:
        return Y

def _plot_barycenter(X, Y, matchings):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # n_y = len(Y.points)
    for i in range(len(X)):
        indices = matchings[i]
        n_i = len(X[i])

        for (y_j, x_i_j) in enumerate(indices):
            y = Y[y_j]
            if y[0] != y[1]:
                if x_i_j < n_i:  # not mapped with the diag
                    x = X[i][x_i_j]
                else:  # y_j is matched to the diagonal
                    x = _proj_on_diag(y)
                ax.plot([y[0], x[0]], [y[1], x[1]], c='black',
                        linestyle="dashed")
    
    ax.scatter(Y[:,0], Y[:,1], color='purple', marker='d')

    for dgm in X:
        ax.scatter(dgm[:,0], dgm[:,1], marker ='o')

    shift = 0.1  # for improved rendering
    xmin = min([np.min(x[:,0]) for x in X]) - shift
    xmax = max([np.max(x[:,0]) for x in X]) + shift
    ymin = min([np.max(x[:,1]) for x in X]) - shift
    ymax = max([np.max(x[:,1]) for x in X]) + shift
    themin = min(xmin, ymin)
    themax = max(xmax, ymax)
    ax.set_xlim(themin, themax)
    ax.set_ylim(themin, themax)
    ax.add_patch(Polygon([[themin,themin], [themax,themin], [themax,themax]], fill=True, color='lightgrey'))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal', adjustable='box')
    ax.set_title("example of (estimated) barycenter")


if __name__=="__main__":
    dg1 = np.array([[0.1, 0.12], [0.21, 0.7], [0.4, 0.5], [0.3, 0.4], [0.35, 0.7], [0.5, 0.55], [0.32, 0.42], [0.1, 0.4], [0.2, 0.4]])
    dg2 = np.array([[0.09, 0.11], [0.3, 0.43], [0.5, 0.61], [0.3, 0.7], [0.42, 0.5], [0.35, 0.41], [0.74, 0.9], [0.5, 0.95], [0.35, 0.45], [0.13, 0.48], [0.32, 0.45]])
    dg3 = np.array([[0.1, 0.15], [0.1, 0.7], [0.2, 0.22], [0.55, 0.84], [0.11, 0.91], [0.61, 0.75], [0.33, 0.46], [0.12, 0.41], [0.32, 0.48]])
    X = [dg1, dg2, dg3]
    Y, a, e = lagrangian_barycenter(X, verbose=True)
