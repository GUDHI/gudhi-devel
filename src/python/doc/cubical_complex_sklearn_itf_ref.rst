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
a TDA scikit-learn pipeline is constructed and is composed of:

#. :class:`~gudhi.sklearn.cubical_persistence.CubicalPersistence` that builds a cubical complex from the inputs and
   returns its persistence diagrams
#. :class:`~gudhi.representations.preprocessing.DiagramSelector` that removes non-finite persistence diagrams values
#. :class:`~gudhi.representations.vector_methods.PersistenceImage` that builds the persistence images from persistence diagrams
#. `SVC <https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html>`_ which is a scikit-learn support
   vector classifier.

This ML pipeline is trained to detect if the hand written digit is an '8' or not, thanks to the fact that an '8' has
two holes in :math:`\mathbf{H}_1`, or, like in this example, three connected components in :math:`\mathbf{H}_0`.

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
    
    X, y = fetch_openml("mnist_784", version=1, return_X_y=True, as_frame=False)
    
    # Target is: "is an eight ?"
    y = (y == "8") * 1
    print("There are", np.sum(y), "eights out of", len(y), "numbers.")
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=0)
    pipe = Pipeline(
        [
            ("cub_pers", CubicalPersistence(homology_dimensions=0, newshape=[-1, 28, 28], n_jobs=-2)),
            # Or for multiple persistence dimension computation
            # ("cub_pers", CubicalPersistence(homology_dimensions=[0, 1], newshape=[-1, 28, 28])),
            # ("H0_diags", DimensionSelector(index=0), # where index is the index in homology_dimensions array
            ("finite_diags", DiagramSelector(use=True, point_type="finite")),
            (
                "pers_img",
                PersistenceImage(bandwidth=50, weight=lambda x: x[1] ** 2, im_range=[0, 256, 0, 256], resolution=[20, 20]),
            ),
            ("svc", SVC()),
        ]
    )
    
    # Learn from the train subset
    pipe.fit(X_train, y_train)
    # Predict from the test subset
    predicted = pipe.predict(X_test)
    
    print(f"Classification report for TDA pipeline {pipe}:\n" f"{metrics.classification_report(y_test, predicted)}\n")

.. code-block:: none

    There are 6825 eights out of 70000 numbers.
    Classification report for TDA pipeline Pipeline(steps=[('cub_pers',
                     CubicalPersistence(newshape=[28, 28], n_jobs=-2)),
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

Cubical complex persistence scikit-learn like interface reference
-----------------------------------------------------------------

.. autoclass:: gudhi.sklearn.cubical_persistence.CubicalPersistence
   :members:
   :special-members: __init__
   :show-inheritance: