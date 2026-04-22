:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Steenrod persistence user manual
=================================

.. include:: steenrod_persistence_sum.inc

Steenrod squares
----------------

Steenrod squares are cohomology operations

.. math::

    \mathrm{Sq}^k : H^n(X;\mathbb{F}_2) \longrightarrow H^{n+k}(X;\mathbb{F}_2)

defined for all integers :math:`k \geq 0` and all topological spaces :math:`X`.
They are natural transformations — commuting with continuous maps — and together they
determine the cohomology ring structure in a finer sense than the cup product alone.
In particular :math:`\mathrm{Sq}^0 = \mathrm{id}`, and for a class :math:`\alpha \in H^n`
we have :math:`\mathrm{Sq}^n(\alpha) = \alpha \smile \alpha` (the cup-square).
When :math:`k > n`, :math:`\mathrm{Sq}^k` is identically zero.

Persistence Steenrod barcodes
------------------------------

Given a filtered simplicial complex :math:`K_0 \subseteq K_1 \subseteq \cdots \subseteq K_m`
over :math:`\mathbb{F}_2`, ordinary persistent cohomology tracks the birth and death of
cohomology classes across the filtration, yielding a persistence barcode.

The *Sq^k Steenrod barcode* enriches this picture: for each cohomological degree :math:`n`,
it records which persistent cohomology classes in degree :math:`n` have a non-trivial image
under :math:`\mathrm{Sq}^k` in degree :math:`n+k`.  Formally, if a cohomology class
:math:`\alpha` is born at index :math:`b` and dies at index :math:`d`, and the cup-i
cochain representative of :math:`\mathrm{Sq}^k(\alpha)` is non-zero, then the bar
:math:`(b, d)` appears in ``steenrod[n]``.

Two filtered complexes may share identical ordinary barcodes yet differ in their Steenrod
barcodes, making the latter a strictly finer invariant.

Algorithm
---------

The computation proceeds in four stages, following :cite:`lupo2022persistencesteenrod` and
using the fast cup-i product formulas of :cite:`medina2023cupifast`:

1. **Filtration sort** — simplices are sorted by dimension within each filtration value,
   producing a boundary matrix in the correct reduction order.
2. **Persistence reduction** — a standard cohomology twist-reduction (with the clearing
   optimization) reduces the coboundary matrix and records cocycle representatives
   in a triangular matrix :math:`V`.
3. **Steenrod matrix** — for each cocycle representative, the cup-i product formula is
   applied to produce a cochain representing :math:`\mathrm{Sq}^k` of that class.
   This stage is parallelised with OpenMP.
4. **Steenrod reduction** — the augmented coboundary-plus-Steenrod matrix is reduced to
   extract the Steenrod barcode.  This stage is also parallelised with OpenMP.

No prior call to :func:`~gudhi.SimplexTree.compute_persistence` is required.

Python interface
-----------------

.. autoclass:: gudhi.SimplexTree
   :members: compute_steenrod_barcodes
   :noindex:

.. autoclass:: gudhi.simplex_tree.SteenrodBarcodes
   :members:
   :noindex:

Example — :math:`\mathbb{R}P^2`
---------------------------------

The real projective plane :math:`\mathbb{R}P^2` is a classical test case:
:math:`H^*(\mathbb{R}P^2;\mathbb{F}_2) \cong \mathbb{F}_2` in each degree 0, 1, and 2,
and :math:`\mathrm{Sq}^1 : H^1 \to H^2` is an isomorphism.  Therefore the Sq\ :sup:`1`
Steenrod barcode contains exactly one essential bar, in degree 2 (tracking the class
detected via Sq\ :sup:`1` applied to the degree-1 generator).

.. testcode::

    import math
    from gudhi import SimplexTree

    # Minimal triangulation of RP^2: 6 vertices, 15 edges, 10 triangles.
    rp2_top = [
        [1, 2, 4], [2, 3, 4], [1, 3, 5], [2, 3, 5], [1, 4, 5],
        [1, 2, 6], [1, 3, 6], [3, 4, 6], [2, 5, 6], [4, 5, 6],
    ]

    st = SimplexTree()
    for tri in rp2_top:
        st.insert(tri, filtration=0.0)

    bars = st.compute_steenrod_barcodes(k=1)

    # Ordinary barcode: H0 and H1 each have one essential bar.
    assert len([b for b in bars.ordinary[0] if math.isinf(b[1])]) == 1
    assert len([b for b in bars.ordinary[1] if math.isinf(b[1])]) == 1

    # Steenrod barcode: exactly one essential bar in degree 2.
    essential_sq1 = [b for b in bars.steenrod[2] if math.isinf(b[1])]
    print("Essential Sq^1 bars in degree 2:", len(essential_sq1))

.. testoutput::

    Essential Sq^1 bars in degree 2: 1
