.. list-table::
   :width: 100%
   :header-rows: 0

   * - :Since: GUDHI 3.6.0
     - :License: MIT
     - :Requires: `TensorFlow <installation.html#tensorflow>`_

We provide TensorFlow 2 models that can handle automatic differentiation for the computation of persistence diagrams from complexes available in the Gudhi library. 
This includes simplex trees, cubical complexes and Vietoris-Rips complexes. Detailed example on how to use these layers in practice are available
in the following `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-optimization.ipynb>`_. Note that even if TensorFlow GPU is enabled, all 
internal computations using Gudhi will be done on CPU. 
