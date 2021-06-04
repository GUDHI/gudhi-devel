:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Cubical complex user manual
===========================
Definition
----------

.. list-table::
   :widths: 25 50 25
   :header-rows: 0

   * - :Author: Pawel Dlotko
     - :Since: GUDHI 2.0.0
     - :License: MIT
   * - :doc:`cubical_complex_user`
     - * :doc:`cubical_complex_ref`
       * :doc:`periodic_cubical_complex_ref`
       * :doc:`cubical_complex_sklearn_itf_ref`
     -


The cubical complex is an example of a structured complex useful in computational mathematics (specially rigorous
numerics) and image analysis.

An *elementary interval* is an interval of a form :math:`[n,n+1]`, or :math:`[n,n]`, for :math:`n \in \mathcal{Z}`.
The first one is called *non-degenerate*, while the second one is a *degenerate* interval. A
*boundary of a elementary interval* is a chain  :math:`\partial [n,n+1] = [n+1,n+1]-[n,n]` in case of
non-degenerated elementary interval and :math:`\partial [n,n] = 0` in case of degenerate elementary interval. An
*elementary cube* :math:`C` is a product of elementary intervals, :math:`C=I_1 \times \ldots \times I_n`.
*Embedding dimension* of a cube is n, the number of elementary intervals (degenerate or not) in the product.
A *dimension of a cube* :math:`C=I_1 \times ... \times I_n` is the number of non degenerate elementary
intervals in the product. A *boundary of a cube* :math:`C=I_1 \times \ldots \times I_n` is a chain obtained
in the following way:

.. math::

    \partial C = (\partial I_1 \times \ldots \times I_n) + (I_1 \times \partial I_2 \times \ldots \times I_n) +
    \ldots + (I_1 \times I_2 \times \ldots \times \partial I_n).

A *cubical complex* :math:`\mathcal{K}` is a collection of cubes closed under operation of taking boundary
(i.e. boundary of every cube from the collection is in the collection). A cube :math:`C` in cubical complex
:math:`\mathcal{K}` is *maximal* if it is not in a boundary of any other cube in :math:`\mathcal{K}`. A
*support* of a cube :math:`C` is the set in :math:`\mathbb{R}^n` occupied by :math:`C` (:math:`n` is the embedding
dimension of :math:`C`).

Cubes may be equipped with a filtration values in which case we have filtered cubical complex. All the cubical
complexes considered in this implementation are filtered cubical complexes (although, the range of a filtration may
be a set of two elements).

For further details and theory of cubical complexes, please consult :cite:`kaczynski2004computational` as well as the
following paper :cite:`peikert2012topological`.

Data structure
--------------

The implementation of Cubical complex provides a representation of complexes that occupy a rectangular region in
:math:`\mathbb{R}^n`. This extra assumption allows for a memory efficient way of storing cubical complexes in a form
of so called bitmaps. Let
:math:`R = [b_1,e_1] \times \ldots \times [b_n,e_n]`, for :math:`b_1,...b_n,e_1,...,e_n \in \mathbb{Z}`,
:math:`b_i \leq d_i` be the considered rectangular region and let :math:`\mathcal{K}` be a filtered
cubical complex having the rectangle :math:`R` as its support. Note that the structure of the coordinate system gives
a way a lexicographical ordering of cells of :math:`\mathcal{K}`. This ordering is a base of the presented
bitmap-based implementation. In this implementation, the whole cubical complex is stored as a vector of the values
of filtration. This, together with dimension of :math:`\mathcal{K}` and the sizes of :math:`\mathcal{K}` in all
directions, allows to determine, dimension, neighborhood, boundary and coboundary of every cube
:math:`C \in \mathcal{K}`.

.. figure::
    ../../doc/Bitmap_cubical_complex/Cubical_complex_representation.png
    :alt: Cubical complex.
    :figclass: align-center

    Cubical complex.

Note that the cubical complex in the figure above is, in a natural way, a product of one dimensional cubical
complexes in :math:`\mathbb{R}`. The number of all cubes in each direction is equal :math:`2n+1`, where :math:`n` is
the number of maximal cubes in the considered direction. Let us consider a cube at the position :math:`k` in the
bitmap.
Knowing the sizes of the bitmap, by a series of modulo operation, we can determine which elementary intervals are
present in the product that gives the cube :math:`C`. In a similar way, we can compute boundary and the coboundary of
each cube. Further details can be found in the literature.

Input Format
------------

In the current implantation, filtration is given at the maximal cubes, and it is then extended by the lower star
filtration to all cubes. There are a number of constructors that can be used to construct cubical complex by users
who want to use the code directly. They can be found in the :doc:`cubical_complex_ref`.
Currently one input from a text file is used. It uses a format inspired from the Perseus software
`Perseus software <http://www.sas.upenn.edu/~vnanda/perseus/>`_ by Vidit Nanda.

.. note::
    While Perseus assume the filtration of all maximal cubes to be non-negative, over here we do not enforce this and
    we allow any filtration values. As a consequence one cannot use ``-1``'s to indicate missing cubes. If you have
    missing cubes in your complex, please set their filtration to :math:`+\infty` (aka. ``inf`` in the file).

The file format is described in details in `Perseus file format <fileformats.html#perseus>`_ section.

.. testcode::

    import gudhi
    cubical_complex = gudhi.CubicalComplex(perseus_file=gudhi.__root_source_dir__ + \
        '/data/bitmap/cubicalcomplexdoc.txt')
    result_str = 'Cubical complex is of dimension ' + repr(cubical_complex.dimension()) + ' - ' + \
        repr(cubical_complex.num_simplices()) + ' simplices.'
    print(result_str)

the program output is:

.. testoutput::
    
    Cubical complex is of dimension 2 - 49 simplices.

Periodic boundary conditions
----------------------------

Often one would like to impose periodic boundary conditions to the cubical complex (cf.
:doc:`periodic_cubical_complex_ref`).
Let :math:`I_1\times ... \times I_n` be a box that is decomposed with a cubical complex :math:`\mathcal{K}`.
Imposing periodic boundary conditions in the direction i, means that the left and the right side of a complex
:math:`\mathcal{K}` are considered the same. In particular, if for a bitmap :math:`\mathcal{K}` periodic boundary
conditions are imposed in all directions, then complex :math:`\mathcal{K}` became n-dimensional torus. One can use
various constructors from the file Bitmap_cubical_complex_periodic_boundary_conditions_base.h to construct cubical
complex with periodic boundary conditions.

One can also use Perseus style input files (see `Perseus file format <fileformats.html#perseus>`_) for the specific periodic case:

.. testcode::

    import gudhi
    periodic_cc = gudhi.PeriodicCubicalComplex(perseus_file=gudhi.__root_source_dir__ + \
        '/data/bitmap/periodiccubicalcomplexdoc.txt')
    result_str = 'Periodic cubical complex is of dimension ' + repr(periodic_cc.dimension()) + ' - ' + \
        repr(periodic_cc.num_simplices()) + ' simplices.'
    print(result_str)

the program output is:

.. testoutput::
    
    Periodic cubical complex is of dimension 2 - 42 simplices.

Or it can be defined as follows:

.. testcode::

    from gudhi import PeriodicCubicalComplex as pcc
    periodic_cc = pcc(top_dimensional_cells = [[0, 0, 0], [0, 1, 0], [0, 0, 0]],
         periodic_dimensions=[True, False])
    result_str = 'Periodic cubical complex is of dimension ' + repr(periodic_cc.dimension()) + ' - ' + \
        repr(periodic_cc.num_simplices()) + ' simplices.'
    print(result_str)

the program output is:

.. testoutput::

    Periodic cubical complex is of dimension 2 - 42 simplices.

Examples
--------

End user programs are available in python/example/ folder.

Tutorial
--------

This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-cubical-complexes.ipynb>`_
explains how to represent sublevels sets of functions using cubical complexes.

Scikit-learn like interface example
-----------------------------------

.. code-block:: python

    # Standard scientific Python imports
    import numpy as np
    # Standard scikit-learn imports
    from sklearn.datasets import fetch_openml
    from sklearn.pipeline import Pipeline
    from sklearn.model_selection import train_test_split
    from sklearn.svm import SVC
    from sklearn import metrics
    
    # Import TDA pipeline requirements
    from gudhi.sklearn.cubical_persistence import CubicalPersistence
    from gudhi.representations import PersistenceImage, DiagramSelector
    
    X, y = fetch_openml('mnist_784', version=1, return_X_y=True, as_frame=False)
    
    # Target is: "is an eight ?"
    y = (y == '8') * 1
    print('There are', np.sum(y), 'eights out of', len(y), 'numbers.')
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=0)
    pipe = Pipeline([('cub_pers', CubicalPersistence(persistence_dim = 0, dimensions=[28,28], n_jobs=-2)),
                     ('finite_diags', DiagramSelector(use=True, point_type="finite")),
                     ('pers_img', PersistenceImage(bandwidth=50,
                                                   weight=lambda x: x[1]**2,
                                                   im_range=[0,256,0,256],
                                                   resolution=[20,20])),
                     ('svc', SVC())])
    
    predicted = pipe.predict(X_test)
    
    print(f"Classification report for TDA pipeline {pipe}:\n"
          f"{metrics.classification_report(y_test, predicted)}\n")

.. code-block:: none

    There are 6825 eights out of 70000 numbers.
    Classification report for TDA pipeline Pipeline(steps=[('cub_pers',
                     CubicalPersistence(dimensions=[28, 28], n_jobs=-2)),
                    ('finite_diags', DiagramSelector(use=True)),
                    ('pers_img',
                     PersistenceImage(bandwidth=50, im_range=[0, 256, 0, 256],
                                      weight=<function <lambda> at 0x7f3e54137ae8>)),
                    ('svc', SVC())]):
                  precision    recall  f1-score   support
    
               0       0.97      0.99      0.98     25284
               1       0.92      0.68      0.78      2716
    
        accuracy                           0.96     28000
       macro avg       0.94      0.84      0.88     28000
    weighted avg       0.96      0.96      0.96     28000