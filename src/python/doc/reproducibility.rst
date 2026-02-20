:orphan:

.. To get rid of WARNING: document isn't included in any toctree

===============
Reproducibility
===============

Some of the GUDHI functionnalities are using randomness. In order to reproduce the results, a specific GUDHI random
generator is available.

This generator inherits from
`numpy random BitGenerator <https://numpy.org/doc/stable/reference/random/bit_generators/generated/numpy.random.BitGenerator.html>`_
and can be constructed with (for reproducibility) or without (for randomness) a seed.

This generator needs to also be able to set CGAL randomness, and this is the reason why the GUDHI random generator is
using internally the
`CGAL random generator <https://doc.cgal.org/latest/Generator/classCGAL_1_1Random.html>`_.

GudhiBitGenerator
=================

.. table::
   :widths: 50 50

   +---------------------------------------------+------------------------------------------+
   | :Requires: `CGAL <installation.html#cgal>`_ | :License: MIT (`LGPL v3 </licensing/>`_) |
   +---------------------------------------------+------------------------------------------+

.. autoclass:: gudhi.random.GudhiBitGenerator
    :members:
    :show-inheritance:

Here is an example of how to use it:

.. code-block:: python

    from numpy.random import Generator
    from gudhi.random import GudhiBitGenerator
    
    rng = Generator(GudhiBitGenerator(42)) # Sets the seed
    rng.random(5)
    # Outputs: array([0.34270148, 0.11108528, 0.42233896, 0.08111117, 0.85644071])


.. note::
    
    This is a work in progress, and for the moment, only `datasets generators <datasets.html#datasets-generators>`_
    can accept a `GudhiBitGenerator` as an argument.
