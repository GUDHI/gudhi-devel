:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Delaunay complex reference manual
*********************************

.. autoclass:: gudhi.filtrations.DelaunayComplex
   :members:
   :inherited-members:
   :exclude-members: __new__, get_point, get_float_relative_precision, set_float_relative_precision

Static methods
==============

.. static nanobind methods fails to be rendered automatically with Sphinx
 cf. https://github.com/wjakob/nanobind/discussions/707
 A workaround is to exclude them from autoclass and to use automethod on it

.. automethod:: gudhi.filtrations.DelaunayComplex.get_float_relative_precision
.. automethod:: gudhi.filtrations.DelaunayComplex.set_float_relative_precision

Alpha complex reference manual
******************************

.. autoclass:: gudhi.filtrations.AlphaComplex
   :members:
   :inherited-members:
   :exclude-members: __new__, get_float_relative_precision, set_float_relative_precision
   :show-inheritance:

Delaunay Čech complex reference manual
**************************************

.. autoclass:: gudhi.filtrations.DelaunayCechComplex
   :members:
   :inherited-members:
   :exclude-members: __new__, get_point, get_float_relative_precision, set_float_relative_precision
   :show-inheritance:
