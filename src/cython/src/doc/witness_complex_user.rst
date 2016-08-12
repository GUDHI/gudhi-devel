===========================
Witness complex user manual
===========================
Definition
----------

.. include:: witness_complex_sum.rst

Implementation
--------------

The principal class of this module is Gudhi::Witness_complex.

In both cases, the constructor for this class takes a {witness}x{closest_landmarks} table, where each row represents a
witness and consists of landmarks sorted by distance to this witness.

.. todo::
    This table can be constructed by two additional classes Landmark_choice_by_furthest_point and
    Landmark_choice_by_random_point also included in the module.

.. figure::
    img/bench_Cy8.png
    :align: center
    
    Running time as function on number of landmarks.

.. figure::
    img/bench_sphere.png
    :align: center
    
    Running time as function on number of witnesses for |L|=300.
