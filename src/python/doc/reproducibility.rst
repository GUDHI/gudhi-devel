:orphan:

.. To get rid of WARNING: document isn't included in any toctree

===============
Reproducibility
===============

Some of the GUDHI functionnalities are using randomness. In order to reproduce the results, a seed mechanism is
available.

.. autofunction:: gudhi.random.set_seed

.. note::

    :func:`~gudhi.random.set_seed` does not set NumPy seed, while some GUDHI functionnalities are using NumPy random
    methods. `NumPy seed <https://numpy.org/doc/stable/reference/random/generated/numpy.random.seed.html>`_ can be set
    with:
    
    .. code-block:: python
    
       import numpy as np       
       np.random.seed(42)
