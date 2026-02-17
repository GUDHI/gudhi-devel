:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Vineyard user manual
====================
.. include:: vineyard_sum.inc

Definition
----------

TODO


Simple example with general filtration
--------------------------------------

This example builds the vineyard from a hand constructed filtration.

.. plot::
    :context:
    :include-source:

    import numpy as np
    from gudhi.vineyard import Vineyard
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

