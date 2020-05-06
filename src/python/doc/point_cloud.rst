:orphan:

.. To get rid of WARNING: document isn't included in any toctree

============================
Point cloud utilities manual
============================

File Readers
------------

.. autofunction:: gudhi.read_points_from_off_file

.. autofunction:: gudhi.read_lower_triangular_matrix_from_csv_file

Subsampling
-----------

:Requires: :ref:`Eigen` :math:`\geq` 3.1.0 and :ref:`CGAL` :math:`\geq` 4.11.0

.. automodule:: gudhi.subsampling
   :members:
   :special-members:
   :show-inheritance:

Time Delay Embedding
--------------------

.. autoclass:: gudhi.point_cloud.timedelay.TimeDelayEmbedding
   :members:
   :special-members: __call__

K nearest neighbors
-------------------

:Requires: :ref:`PyKeOps`, :ref:`SciPy`, :ref:`Scikit-learn`, and/or :ref:`Hnswlib` in function of the selected `implementation`.

.. automodule:: gudhi.point_cloud.knn
   :members:
   :undoc-members:
   :special-members: __init__

Distance to measure
-------------------

.. automodule:: gudhi.point_cloud.dtm
   :members:
   :undoc-members:
   :special-members: __init__
