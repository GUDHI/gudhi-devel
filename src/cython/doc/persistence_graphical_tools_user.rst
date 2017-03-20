=======================================
Persistence graphical tools user manual
=======================================
Definition
----------
.. include:: persistence_graphical_tools_sum.rst


Show palette values
-------------------

This function is useful to show the color palette values of dimension:


.. testcode::

    import gudhi
    gudhi.show_palette_values(alpha=1.0)

.. plot::

    import gudhi
    gudhi.show_palette_values(alpha=1.0)

Show persistence as a barcode
-----------------------------

This function can display the persistence result as a barcode:

.. testcode::

    import gudhi
    
    periodic_cc = gudhi.PeriodicCubicalComplex(perseus_file='3d_torus.txt')
    diag = periodic_cc.persistence()
    gudhi.plot_barcode_persistence(diag)

.. plot::

    import gudhi

    periodic_cc = gudhi.PeriodicCubicalComplex(perseus_file='3d_torus.txt')
    diag = periodic_cc.persistence()
    gudhi.plot_barcode_persistence(diag)

Show persistence as a diagram
-----------------------------

This function can display the persistence result as a diagram:

.. testcode::

    import gudhi
    
    rips_complex = gudhi.RipsComplex(off_file='tore3D_300.off', max_edge_length=2.0)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=3)
    diag = simplex_tree.persistence()
    gudhi.plot_diagram_persistence(diag)

.. plot::

    import gudhi

    rips_complex = gudhi.RipsComplex(off_file='tore3D_300.off', max_edge_length=2.0)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=3)
    diag = simplex_tree.persistence()
    gudhi.plot_diagram_persistence(diag)
