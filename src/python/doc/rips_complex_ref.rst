:orphan:

.. To get rid of WARNING: document isn't included in any toctree

=============================
Rips complex reference manual
=============================

.. autoclass:: gudhi.RipsComplex
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: gudhi.RipsComplex.__init__

======================================
Weighted Rips complex reference manual
======================================

.. autoclass:: gudhi.weighted_rips_complex.WeightedRipsComplex
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: gudhi.weighted_rips_complex.WeightedRipsComplex.__init__

Basic examples
--------------

The following example computes the weighted Rips filtration associated with a distance matrix and weights on vertices.

.. testcode::

    from gudhi.weighted_rips_complex import WeightedRipsComplex
    dist = [[], [1]]
    weights = [1, 100]
    w_rips = WeightedRipsComplex(distance_matrix=dist, weights=weights)
    st = w_rips.create_simplex_tree(max_dimension=2)
    print(st.get_filtration())

The output is:

.. testoutput::

    [([0], 2.0), ([1], 200.0), ([0, 1], 200.0)]

Combining with DistanceToMeasure, one can compute the DTM-filtration of a point set, as in `this notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-DTM-filtrations.ipynb>`_. 

.. testcode::

    import numpy as np
    from scipy.spatial.distance import cdist
    from gudhi.point_cloud.dtm import DistanceToMeasure
    from gudhi.weighted_rips_complex import WeightedRipsComplex
    pts = np.array([[2.0, 2.0], [0.0, 1.0], [3.0, 4.0]])
    dist = cdist(pts,pts)
    dtm = DistanceToMeasure(2, q=2, metric="precomputed")
    r = dtm.fit_transform(dist)
    w_rips = WeightedRipsComplex(distance_matrix=dist, weights=r)
    st = w_rips.create_simplex_tree(max_dimension=2)
    print(st.persistence())

The output is:

.. testoutput::

    [(0, (3.1622776601683795, inf)), (0, (3.1622776601683795, 5.39834563766817)), (0, (3.1622776601683795, 5.39834563766817))]
