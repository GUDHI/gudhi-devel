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

:Requires: `Eigen <installation.html#eigen>`_ :math:`\geq` 3.1.0 and `CGAL <installation.html#cgal>`_ :math:`\geq` 4.11.0

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
