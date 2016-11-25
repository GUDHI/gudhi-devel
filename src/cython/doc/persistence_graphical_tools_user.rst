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

Show persistence as a barcode
-----------------------------

This function can display the persistence result as a barcode:

.. testcode::

    import gudhi
    
    periodic_cc = gudhi.PeriodicCubicalComplex(perseus_file='3d_torus.txt')
    diag = periodic_cc.persistence()
    gudhi.barcode_persistence(diag)


Show persistence as a diagram
-----------------------------

This function can display the persistence result as a diagram:

.. testcode::

    import gudhi
    
    alpha_complex = gudhi.AlphaComplex(off_file='tore3D_300.off')
    simplex_tree = gudhi.SimplexTree()
    alpha_complex.create_simplex_tree(simplex_tree)
    diag = simplex_tree.persistence()
    gudhi.diagram_persistence(diag)
