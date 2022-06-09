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

We provide scikit-learn classes for computing cover complexes, i.e., Nerve, Graph Induced and Mapper complexes. Detailed examples on how to use these classes in practice are available
in the following `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-cover-complex.ipynb>`_. 

Example of Mapper cover complex computed from a point cloud
-----------------------------------------------------------

.. testcode::

    import numpy as np
    from sklearn.cluster import AgglomerativeClustering
    from gudhi.sklearn import MapperComplex

    X = np.array([[1,1],[1,1.5],[1,2],[1,2.5],[1,3],[1.5,2],[1.5,3],[2,2],[2,2.5],[2,3],[2,3.5],[2,4]])
    F = np.array([[1,1,1,1,1,1.5,1.5,2,2,2,2,2],[1,1.5,2,2.5,3,2,3,2,2.5,3,3.5,4]]).T


    Mapper = MapperComplex(
        input_type="point cloud",
        colors=None,
        mask=0,
        filters=F,
        filter_bnds=np.array([[0.5, 2.5], [0.5, 4.5]]),
        resolutions=np.array([2, 4]),
        gains=np.array([0.3, 0.3]),
        clustering=AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=0.6),
    )
    
    Mapper.fit(X)

    print([s for s,_ in Mapper.simplex_tree.get_simplices()])

.. testoutput::

    [[0, 1], [0], [1, 2], [1, 3], [1], [2, 4], [2], [3, 4], [3], [4, 5], [4], [5]]
    
Documentation for MapperComplex
-------------------------------

.. autoclass:: gudhi.sklearn.MapperComplex
   :members:
   :special-members: __init__
   :show-inheritance:

Documentation for GraphInducedComplex
-------------------------------------

.. autoclass:: gudhi.sklearn.GraphInducedComplex
   :members:
   :special-members: __init__
   :show-inheritance:

Documentation for NerveComplex
------------------------------

.. autoclass:: gudhi.sklearn.NerveComplex
   :members:
   :special-members: __init__
   :show-inheritance:
