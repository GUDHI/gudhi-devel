Persistent cohomology user manual
=================================
Definition
----------
=====================================  =====================================  =====================================
:Author: Clément Maria                 :Introduced in: GUDHI PYTHON 2.0.0     :Copyright: GPL v3
=====================================  =====================================  =====================================

+---------------------------------------------+----------------------------------------------------------------------+
|  :doc:`persistent_cohomology_user`          | Please refer to each data structure that contains persistence        |
|                                             | feature for reference:                                               |
|                                             |                                                                      |
|                                             | * :doc:`simplex_tree_ref`                                            |
+---------------------------------------------+----------------------------------------------------------------------+


Computation of persistent cohomology using the algorithm of :cite:`DBLP:journals/dcg/SilvaMV11` and
:cite:`DBLP:journals/corr/abs-1208-5018` and the Compressed Annotation Matrix implementation of
:cite:`DBLP:conf/esa/BoissonnatDM13`.
     
The theory of homology consists in attaching to a topological space a sequence of (homology) groups, capturing global
topological features like connected components, holes, cavities, etc. Persistent homology studies the evolution --
birth, life and death -- of these features when the topological space is changing. Consequently, the theory is
essentially composed of three elements:

* topological spaces
* their homology groups
* an evolution scheme.

Topological Spaces
------------------

Topological spaces are represented by simplicial complexes.
Let :math:`V = \{1, \cdots ,|V|\}` be a set of *vertices*.
A *simplex* :math:`\sigma` is a subset of vertices :math:`\sigma \subseteq V`.
A *simplicial complex* :math:`\mathbf{K}` on :math:`V` is a collection of simplices :math:`\{\sigma\}`,
:math:`\sigma \subseteq V`, such that :math:`\tau \subseteq \sigma \in \mathbf{K} \Rightarrow \tau \in \mathbf{K}`.
The dimension :math:`n=|\sigma|-1` of :math:`\sigma` is its number of elements minus 1.
A *filtration* of a simplicial complex is a function :math:`f:\mathbf{K} \rightarrow \mathbb{R}` satisfying
:math:`f(\tau)\leq f(\sigma)` whenever :math:`\tau \subseteq \sigma`.

Homology
--------

For a ring :math:`\mathcal{R}`, the group of *n-chains*, denoted :math:`\mathbf{C}_n(\mathbf{K},\mathcal{R})`, of
:math:`\mathbf{K}` is the group of formal sums of n-simplices with :math:`\mathcal{R}` coefficients. The
*boundary operator* is a linear operator
:math:`\partial_n: \mathbf{C}_n(\mathbf{K},\mathcal{R}) \rightarrow \mathbf{C}_{n-1}(\mathbf{K},\mathcal{R})`
such that :math:`\partial_n \sigma = \partial_n [v_0, \cdots , v_n] = \sum_{i=0}^n (-1)^{i}[v_0,\cdots ,\widehat{v_i}, \cdots,v_n]`,
where :math:`\widehat{v_i}` means :math:`v_i` is omitted from the list. The chain groups form a sequence:

.. math::

    \cdots \ \ \mathbf{C}_n(\mathbf{K},\mathcal{R}) \xrightarrow{\ \partial_n\ }
    \mathbf{C}_{n-1}(\mathbf{K},\mathcal{R}) \xrightarrow{\partial_{n-1}} \cdots \xrightarrow{\ \partial_2 \ }
    \mathbf{C}_1(\mathbf{K},\mathcal{R}) \xrightarrow{\ \partial_1 \ }  \mathbf{C}_0(\mathbf{K},\mathcal{R})

of finitely many groups :math:`\mathbf{C}_n(\mathbf{K},\mathcal{R})` and homomorphisms :math:`\partial_n`, indexed by
the dimension :math:`n \geq 0`. The boundary operators satisfy the property :math:`\partial_n \circ \partial_{n+1}=0`
for every :math:`n > 0` and we define the homology groups:

.. math::

    \mathbf{H}_n(\mathbf{K},\mathcal{R}) = \ker \partial_n / \mathrm{im} \  \partial_{n+1}

We refer to :cite:`Munkres-elementsalgtop1984` for an introduction to homology
theory and to :cite:`DBLP:books/daglib/0025666` for an introduction to persistent homology.

Indexing Scheme
---------------

"Changing" a simplicial complex consists in applying a simplicial map. An *indexing scheme* is a directed graph
together with a traversal order, such that two consecutive nodes in the graph are connected by an arrow (either forward
or backward).
The nodes represent simplicial complexes and the directed edges simplicial maps.

From the computational point of view, there are two types of indexing schemes of interest in persistent homology:
    
* linear ones
  :math:`\bullet \longrightarrow \bullet \longrightarrow \cdots \longrightarrow \bullet \longrightarrow \bullet`
  in persistent homology :cite:`DBLP:journals/dcg/ZomorodianC05`,
* zigzag ones
  :math:`\bullet \longrightarrow \bullet \longleftarrow \cdots \longrightarrow \bullet \longleftarrow \bullet`
  in zigzag persistent homology :cite:`DBLP:journals/focm/CarlssonS10`.
  
These indexing schemes have a natural left-to-right traversal order, and we describe them with ranges and iterators.
In the current release of the Gudhi library, only the linear case is implemented.

In the following, we consider the case where the indexing scheme is induced by a filtration.

Ordering the simplices by increasing filtration values (breaking ties so as a simplex appears after its subsimplices of
same filtration value) provides an indexing scheme.

Examples
--------

We provide several example files: run these examples with -h for details on their use.

.. only:: builder_html

    * :download:`alpha_complex_diagram_persistence_from_off_file_example.py <../example/alpha_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`periodic_cubical_complex_barcode_persistence_from_perseus_file_example.py <../example/periodic_cubical_complex_barcode_persistence_from_perseus_file_example.py>`
    * :download:`rips_complex_diagram_persistence_from_off_file_example.py <../example/rips_complex_diagram_persistence_from_off_file_example.py>`
    * :download:`rips_persistence_diagram.py <../example/rips_persistence_diagram.py>`
    * :download:`rips_complex_diagram_persistence_from_distance_matrix_file_example.py <../example/rips_complex_diagram_persistence_from_distance_matrix_file_example.py>`
    * :download:`random_cubical_complex_persistence_example.py <../example/random_cubical_complex_persistence_example.py>`
    * :download:`tangential_complex_plain_homology_from_off_file_example.py <../example/tangential_complex_plain_homology_from_off_file_example.py>`

Bibliography
============

.. bibliography:: bibliography.bib
   :filter: docnames
   :style: unsrt
