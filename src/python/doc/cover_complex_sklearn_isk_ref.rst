:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Scikit-learn class for cover complexes
######################################

.. list-table::
   :widths: 40 30 30
   :header-rows: 0

   * - :Since: GUDHI 3.5.0
     - :License: MIT
     - :Requires: `Scikit-learn <installation.html#scikit-learn>`__

We provide a scikit-learn class for computing cover complexes. Detailed examples on how to use this class in practice are available
in the following `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-cover-complex.ipynb>`_. 

Example of Mapper cover complex computed from a point cloud
-----------------------------------------------------------

.. testcode::

    import numpy as np
    from sklearn.cluster import AgglomerativeClustering
    from gudhi.sklearn import CoverComplex

    X = np.array([[1,1],[1,1.5],[1,2],[1,2.5],[1,3],[1.5,2],[1.5,3],[2,2],[2,2.5],[2,3],[2,3.5],[2,4]])
    F = np.array([[1,1,1,1,1,1.5,1.5,2,2,2,2,2],[1,1.5,2,2.5,3,2,3,2,2.5,3,3.5,4]]).T


    Mapper = CoverComplex(
        complex_type="mapper",
        input_type="point cloud",
        cover="functional",
        colors=None,
        mask=0,
        filters=F,
        filter_bnds=np.array([[0.5, 2.5], [0.5, 4.5]]),
        resolutions=np.array([2, 4]),
        gains=np.array([0.3, 0.3]),
        clustering=AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=0.6),
    )
    
    _ = Mapper.fit(X)

    print(list(Mapper.simplex_tree.get_filtration()))

.. testoutput::

    [([0], -3.0), ([1], -3.0), ([0, 1], -3.0), ([2], -3.0), ([1, 2], -3.0), ([3], -3.0), ([1, 3], -3.0), ([4], -3.0), ([2, 4], -3.0), ([3, 4], -3.0), ([5], -3.0), ([4, 5], -3.0)]
    
Documentation for CoverComplex
------------------------------

.. autoclass:: gudhi.sklearn.CoverComplex
   :members:
   :special-members: __init__
   :show-inheritance:
