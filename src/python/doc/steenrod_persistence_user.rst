:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Steenrod persistence user manual
=================================

.. include:: steenrod_persistence_sum.inc

Steenrod squares
----------------

Steenrod squares are cohomology operations

.. math::

    \mathrm{Sq}^k : H^n(X;\mathbb{F}_2) \longrightarrow
    H^{n+k}(X;\mathbb{F}_2)

defined naturally for all integers :math:`k \geq 0` and all topological
spaces :math:`X`. They provide finer discriminatory information beyond
Betti numbers. For example, with mod 2 coefficients, the real projective
plane and the wedge of a circle and a sphere are indistinguishable using
Betti numbers alone, but :math:`\mathrm{Sq}^k` separates them. With any
coefficients, the complex projective space and the wedge of a 2-sphere
and a 4-sphere are indistinguishable, but :math:`\mathrm{Sq}^k` tells
them apart.

Persistence Steenrod barcodes
-----------------------------

The :math:`\mathrm{Sq}^k` *Steenrod barcode* generalizes the
:math:`\mathbb{F}_2` relative cohomology barcode of a filtered
simplicial complex in the sense that for :math:`k = 0` one recovers it.
It is defined as the barcode of the image of the persistent module
morphism given by
:math:`\mathrm{Sq}^k \colon H^n(X,X_i) \to H^{n+k}(X,X_i)` for each
:math:`i` in the filtration index, and thus records the action of
Steenrod squares along the filtration producing more discriminatory
barcodes. For a more comprehensive treatment of persistent cohomology
operations consult :cite:`medina2025cohomologyoperations`.

Algorithm
---------

The computation proceeds in four stages, following :cite:`lupo2022persistencesteenrod` and
using the fast cup-:math:`i` product formulas of :cite:`medina2023cupifast`:

1. **Filtration sort** — simplices are sorted by dimension within each filtration value,
   producing a boundary matrix in the correct reduction order.
2. **Persistence reduction** — a standard cohomology twist-reduction (with the clearing
   optimization) reduces the coboundary matrix and records cocycle representatives
   in a triangular matrix :math:`V`.
3. **Steenrod matrix** — for each cocycle representative, the cup-:math:`i` product formula is
   applied to produce a cochain representing :math:`\mathrm{Sq}^k` of that class.
   This stage is parallelised with TBB when ``GUDHI_USE_TBB`` is enabled.
4. **Steenrod reduction** — the augmented coboundary-plus-Steenrod matrix is reduced to
   extract the Steenrod barcode.  This stage is also parallelised with TBB.

No prior call to :func:`~gudhi.SimplexTree.compute_persistence` is required.

Python interface
-----------------

.. automethod:: gudhi.SimplexTree.compute_steenrod_barcodes

Example — cone on :math:`\mathbb{C}P^2`
-----------------------------------------

The complex projective plane :math:`\mathbb{C}P^2` satisfies
:math:`H^*(\mathbb{C}P^2;\mathbb{F}_2) \cong \mathbb{F}_2` in degrees 0, 2, and 4,
with :math:`\mathrm{Sq}^2 : H^2 \to H^4` an isomorphism (the generator squares to
the generator).  Coning off :math:`\mathbb{C}P^2` at filtration value 1 kills all
classes, so every bar is finite.  This makes the example a clean illustration
of the relative-cohomology output.

The default **relative** convention (``absolute=False``) returns each
barcode as a pair ``(finites_per_dim, infinites_per_dim)``:

* ``finites[d]`` is a numpy array of shape ``(n_bars, 2)`` whose rows are
  ``(death_value, birth_value)`` with ``death < birth``, encoding the
  relative-cohomology bar ``[a_p, a_{q+1})`` of Lupo, Medina-Mardones,
  Tauzin (2022) §2.4.

* ``infinites[d]`` is a numpy array of shape ``(n_bars,)`` listing the
  birth values of essential bars; the implicit lower endpoint of the
  relative-cohomology interval is ``-inf``.

Sq\ :sup:`k` Steenrod bars sit at one relative dimension higher than the
absolute image degree they represent.

Passing ``absolute=True`` switches to the absolute-cohomology
convention: finite bars at relative dimension ``d`` move down to
absolute dimension ``d - 1`` and essential bars stay where they are.
Numerical values do not change.  The conversion is well-defined only
when the relative ordinary barcode has no essential bars at degrees
``[1, max_dim]`` (the duality bijection of LMT 2022 §3.3); a
``UserWarning`` is emitted otherwise.

.. testcode::

    from gudhi import SimplexTree

    top_cp2 = [
        (1,2,4,5,6),(2,3,5,6,4),(3,1,6,4,5),(1,2,4,5,9),(2,3,5,6,7),
        (3,1,6,4,8),(2,3,6,4,9),(3,1,4,5,7),(1,2,5,6,8),(3,1,5,6,9),
        (1,2,6,4,7),(2,3,4,5,8),(4,5,7,8,9),(5,6,8,9,7),(6,4,9,7,8),
        (4,5,7,8,3),(5,6,8,9,1),(6,4,9,7,2),(5,6,9,7,3),(6,4,7,8,1),
        (4,5,8,9,2),(6,4,8,9,3),(4,5,9,7,1),(5,6,7,8,2),(7,8,1,2,3),
        (8,9,2,3,1),(9,7,3,1,2),(7,8,1,2,6),(8,9,2,3,4),(9,7,3,1,5),
        (8,9,3,1,6),(9,7,1,2,4),(7,8,2,3,5),(9,7,2,3,6),(7,8,3,1,4),
        (8,9,1,2,5),
    ]

    # SimplexTree.insert inserts the simplex and all its faces automatically.
    st = SimplexTree()
    for s in top_cp2:
        st.insert(list(s), filtration=0.0)
    cone_v = 10
    for s in top_cp2:
        st.insert(list(s) + [cone_v], filtration=1.0)

    # Default (relative) convention.  Each of ``ordinary`` and ``steenrod``
    # is a pair ``(finites_per_dim, infinites_per_dim)``.  At k=2 the Sq^2
    # bar lives at relative dimension 5 — one higher than its absolute
    # image degree H^4(CP^2; F_2).
    _, steenrod = st.compute_steenrod_barcodes(k=2)
    st_finites, st_infinites = steenrod
    print("steenrod finites[5]:", st_finites[5].tolist())
    print("steenrod infinites[5]:", st_infinites[5].tolist())

.. testoutput::

    steenrod finites[5]: [[0.0, 1.0]]
    steenrod infinites[5]: []
