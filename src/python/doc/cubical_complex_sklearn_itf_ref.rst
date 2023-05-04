:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Cubical complex persistence scikit-learn like interface
#######################################################

.. list-table::
   :width: 100%
   :header-rows: 0

   * - :Since: GUDHI 3.6.0
     - :License: MIT
     - :Requires: `Scikit-learn <installation.html#scikit-learn>`_

Cubical complex persistence scikit-learn like interface example
---------------------------------------------------------------

In this example, hand written digits are used as an input.
A TDA scikit-learn pipeline is constructed and is composed of:

#. :class:`~gudhi.sklearn.cubical_persistence.CubicalPersistence` that builds a cubical complex from the input images
   and returns its persistence diagrams.
#. :class:`~gudhi.representations.preprocessing.DiagramSelector` that removes non-finite persistence diagrams values.
#. :class:`~gudhi.representations.vector_methods.PersistenceImage` that builds the persistence images from persistence
   diagrams.
#. `SVC <https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html>`_ which is a scikit-learn support
   vector classifier.

This ML pipeline is trained to detect if the hand written digit is an '8' or not, thanks to the fact that an '8' has
two holes in :math:`\mathbf{H}_1`, or, like in this example, three connected components in :math:`\mathbf{H}_0`.

.. literalinclude:: ../../python/example/cubical_complex_sklearn_itf.py
   :language: python

.. code-block:: none

    There are 6825 eights out of 70000 numbers.
    Classification report for TDA pipeline Pipeline(steps=[('cub_pers',
                     CubicalPersistence(homology_dimensions=0, n_jobs=-2)),
                    ('finite_diags', DiagramSelector(use=True)),
                    ('pers_img',
                     PersistenceImage(bandwidth=50, im_range=[0, 256, 0, 256],
                                      weight=<function <lambda> at 0x7fc4578a80d0>)),
                    ('svc', SVC())]):
                  precision    recall  f1-score   support
    
               0       0.97      0.99      0.98     25284
               1       0.92      0.68      0.78      2716
    
        accuracy                           0.96     28000
       macro avg       0.94      0.84      0.88     28000
    weighted avg       0.96      0.96      0.96     28000

Cubical complex persistence scikit-learn like interface reference
-----------------------------------------------------------------

.. autoclass:: gudhi.sklearn.cubical_persistence.CubicalPersistence
   :members:
   :show-inheritance:
