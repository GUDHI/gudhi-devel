:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Vineyard user manual
====================
.. include:: vineyard_sum.inc

Definition
----------

From a filtration of the form
:math:`\mathcal{K}_0 \rightarrow \mathcal{K}_1 \rightarrow \cdots \rightarrow \mathcal{K}_n` we can compute the
persistence diagram. For example, lets take two adjacent triangles and add they faces one after the other:

.. code::

    from gudhi import SimplexTree
    
    st = SimplexTree()
    st.insert([0], 0)
    st.insert([1], 1)
    st.insert([2], 2)
    st.insert([3], 3)
    st.insert([0, 3], 4)
    st.insert([0, 2], 5)
    st.insert([1, 2], 6)
    st.insert([2, 3], 7)
    st.insert([0, 1], 8)
    st.insert([0, 1, 2], 9)
    st.insert([0, 2, 3], 10)

This will result in the following diagram:

.. plot::
    :context: reset

    import numpy as np
    from gudhi import SimplexTree, plot_persistence_diagram
    
    st = SimplexTree()
    st.insert([0], 0)
    st.insert([1], 1)
    st.insert([2], 2)
    st.insert([3], 3)
    st.insert([0, 3], 4)
    st.insert([0, 2], 5)
    st.insert([1, 2], 6)
    st.insert([2, 3], 7)
    st.insert([0, 1], 8)
    st.insert([0, 1, 2], 9)
    st.insert([0, 2, 3], 10)
    plot_persistence_diagram(st.persistence(homology_coeff_field=2, min_persistence=0, persistence_dim_max=True))

But now we could decide to let :math:`[2, 3]` appear before :math:`[1, 2]` and the resulting filtration would still be
valid (any face appears before its cofaces), but the first 1-cycle would be born one step earlier:

.. plot::
    :context: close-figs
    :include-source:
    
    st.assign_filtration([2, 3], 6)
    st.assign_filtration([1, 2], 7)
    plot_persistence_diagram(st.persistence(homology_coeff_field=2, min_persistence=-1, persistence_dim_max=True))

Note that the number of bars does not change as the final complex remains the same. And because persistence diagrams
are stable under the right distance functions, the points in the diagram will not move a lot. Therefore we can "stack"
those diagrams and obtain lines which are tracing the changes of the barcode under those "swap operations".
One of those lines is called a **Vine** and the set of all the vines of the same filtered complex is called a
**Vineyard**.

Implementation Details
----------------------

To ensure the right matching of a point in a persistence diagram with the next diagram and to avoid recomputing
all persistence pairs from scratch, the first layer of the vineyard is computed like any other diagram, but the
remaining layers are computed by updating the columns and rows corresponding to the two swapped cells in the underlying
matrix from the initial computation. The implementation is based on :cite:`vineyards` and :cite:`zigzag`.

The python interface has two classes: :class:`~gudhi.vineyard.Vineyard` and
:class:`~gudhi.vineyard.PointCloudRipsVineyard`. The first one is a general interface, which can take any type of
filtered complex as input. The second is a special case for the first and will take point clouds as input and compute
Rips filtration from them to obtain the vineyards.

The update methods can take filtration values which differs completely from the precedent filtration. The sequence
of swaps will be automatically deduced from them, but the resulting vineyard will only show the start and end state
of the points, not the intermediate swaps.

In addition to the vineyard, the two classes can also retrieve a **representative cycle** for each point in a diagram.

Simple example with general filtration
--------------------------------------

This example builds the vineyard from a hand constructed filtration, based on the same two triangles than above.

.. plot::
    :context: reset
    :include-source:

    import numpy as np
    from gudhi.vineyard import Vineyard

    # boundaries from two adjacent triangles sharing the edge [0, 2]
    # note that different from the simplex tree, the boundary cells are represented by their position in the array
    # and not by their vertices, as the cells do not have to be simplicial
    boundaries = [np.array([]), np.array([]), np.array([]), np.array([]), 
        np.array([0, 3]), np.array([0, 2]), np.array([1, 2]), np.array([2, 3]), np.array([0, 1]), 
        np.array([4, 5, 7]), np.array([5, 6, 8])]
    dimensions = np.array([0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2], dtype=int)
    f1 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=float)
    f2 = np.array([0, 1, 2, 3, 4, 7, 6, 5, 8, 9, 10], dtype=float)
    f3 = np.array([0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 10], dtype=float)

    vy = Vineyard(store_latest_cycles=False)
    vy.initialize(
        boundaries=boundaries, 
        dimensions=dimensions, 
        filtration_values=f1, 
        number_of_updates=2    #not obligatory, but better for memory pre-allocation
    )
    vy.update(filtration_values=f2)
    vy.update(filtration_values=f3)

    vy.plot_vineyards(dim=None, square_scaling=False)


Simple example with rips filtration from file
---------------------------------------------

This example builds the vineyard from rips filtration based on point clouds stored in files.
See :meth:`PointCloudRipsVineyard.from_files <gudhi.vineyard.PointCloudRipsVineyard.from_files>` for more information
about path and file formats.
The data set in this example is provided by
`Francesca Bertoglio <https://github.com/frabertoglio/TDAdynamicdata/tree/main/simulations/rings>`_.

.. plot::
    :context: reset
    :include-source:

    from gudhi.vineyard import PointCloudRipsVineyard

    rvy = PointCloudRipsVineyard.from_files(
        path_prefix="../../data/points/unknotting_rings/rings", 
        store_point_coordinates=True, 
        store_cycles=True
    )

.. list-table::
    :widths: 50 50

    * - .. plot::
                :context: close-figs
                :include-source:

                rvy.plot_vineyards(dim=1, min_bar_length=4.5, noise_option="gray_band")
      - .. plot::
                :context: close-figs
                :include-source:

                rvy.plot_1D_representative_cycles(10, min_bar_length=4.5)

