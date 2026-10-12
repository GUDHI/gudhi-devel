:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Persistence graphical tools user manual
=======================================
Definition
----------
.. include:: persistence_graphical_tools_sum.inc


Show persistence as a barcode
-----------------------------

This function can display the persistence result as a barcode:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import gudhi
    from gudhi.datasets.generators import points
    
    gudhi.random.set_seed(42)
    point_cloud = points.c_2_torus(n_samples=300, major_radius=2., minor_radius=1.)
    
    rips_complex = gudhi.RipsComplex(points=point_cloud, max_edge_length=0.7)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=3)
    diag = simplex_tree.persistence(min_persistence=0.4)
    
    gudhi.plot_persistence_barcode(diag)
    plt.show()


Show persistence as a diagram
-----------------------------

This function can display the persistence result as a diagram:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import numpy as np
    import gudhi
    from gudhi.datasets.generators import points
    
    gudhi.random.set_seed(42)
    point_cloud = points.c_2_torus(n_samples=2500, major_radius=2., minor_radius=1.)
    
    simplex_tree = gudhi.RipsComplex(points=point_cloud, max_edge_length=0.67).create_simplex_tree(max_dimension=3)
    diags = simplex_tree.persistence(min_persistence=0.01)
    # Save the H1 persistence results in a file
    diags_H1 = simplex_tree.persistence_intervals_in_dimension(1)
    np.save("2_torus_H1.npy", diags_H1)
    ax = gudhi.plot_persistence_diagram(diags)
    # We can modify the title, aspect, etc.
    ax.set_title("Persistence diagram of a torus")
    ax.set_aspect("equal")  # forces to be square shaped
    plt.show()


Note that (as barcode and density) it can also take a simple `np.array`
of shape (N x 2) encoding a persistence diagram (in a given dimension).

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import gudhi
    import numpy as np
    d = np.array([[0., 1.], [1., 2.], [1., np.inf]])
    gudhi.plot_persistence_diagram(d)
    plt.show()


Persistence density
-------------------

:Requires: `SciPy <installation.html#scipy>`_

If you want more information on a specific dimension, for instance:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import numpy as np
    import gudhi

    # "2_torus_H1.npy" obtained from np.save method - cf. "Show persistence as a diagram"
    birth_death = np.load("2_torus_H1.npy")
    # Use subplots to display diagram and density side by side
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
    gudhi.plot_persistence_diagram(persistence=birth_death, axes=axes[0])
    gudhi.plot_persistence_density(persistence=birth_death, axes=axes[1])
    plt.show()


LaTeX support
-------------

If you are facing issues with `LaTeX <installation.html#latex>`_ rendering, you can still deactivate LaTeX rendering by
saying:

.. code-block:: python

    import gudhi
    gudhi.persistence_graphical_tools._gudhi_matplotlib_use_tex=False
