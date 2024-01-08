# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Marc Glisse
#
# Copyright (C) 2023 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from .._edge_collapse import _collapse_edges


def reduce_graph(input_edges, nb_iterations=1):
    """Simplify a clique filtration, preserving its persistent homology.
    The clique filtration is represented as a sparse weighted graph, giving for each edge its 2 vertices
    and its filtration value, and the output is in the same format. An edge may be represented arbitrarily as (i, j)
    or (j, i). Listing the same edge twice, or a self-loop, is undefined.
    The cliques of the graph composed of the edges with filtration value less than some arbitrary value implicitly
    define a flag complex, so the weighted graph describes a flag filtration.
    This function outputs another flag filtration which is guaranteed to have the same persistent homology as the input.
    The algorithm works by collapsing edges, as described in :cite:`edgecollapsearxiv`.
    Note that this simplification is independent of the filtration values of the vertices.
    The output has the same shape as the input, which is presumed to be (N, N) where all vertices have index
    less than N, since the simplification does not affect vertices.

    :param input_edges: Input weighted graph.
    :type input_edges: scipy.sparse.coo_matrix
    :param nb_iterations: The number of times we apply the algorithm. Default is 1.
    :type nb_iterations: int
    :rtype: scipy.sparse.coo_matrix
    """
    from scipy.sparse import coo_matrix

    r = _collapse_edges(input_edges.row, input_edges.col, input_edges.data, nb_iterations)
    return coo_matrix((r[1], (r[0][0], r[0][1])), input_edges.shape)
