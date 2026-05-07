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
   This stage is parallelised with OpenMP.
4. **Steenrod reduction** — the augmented coboundary-plus-Steenrod matrix is reduced to
   extract the Steenrod barcode.  This stage is also parallelised with OpenMP.

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
classes, so every bar is finite.  This makes the example ideal for illustrating the
difference between the ``absolute`` and ``relative`` output conventions.

In the **raw relative** convention (``absolute=False``, the default) each bar
is a tuple ``(death_value, birth_value)`` with ``death < birth`` for finite
bars and ``death = -inf`` for essential bars.  This matches the half-open
interval convention of Lupo, Medina-Mardones, Tauzin (2022) §2.4 — the bar
``[p, q]`` represents ``[a_p, a_{q+1})`` with ``a_0 = -inf``.  The slot order
is ``(slot 0 = death, slot 1 = birth)`` and finite bars sit at one degree
higher than the absolute image degree.

In the **absolute** convention (``absolute=True``) bars satisfy
``(birth_value, death_value)`` with ``birth ≤ death``, ``death = +inf`` for
essentials, and finite bars are shifted from relative dimension ``d`` to
absolute dimension ``d - 1``.  This conversion rests on a duality bijection
that requires all ordinary bars at relative dimensions ``[1, max_dim]`` to be
finite; a warning is emitted when this condition is not met.

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

    # Raw relative convention (default): tuples are (death, birth) with
    # death < birth (or death = -inf for essential).  Finite bars sit one
    # degree higher than their absolute image degree.  The Sq^2 bar lives in
    # steenrod[5].
    _, steenrod_rel = st.compute_steenrod_barcodes(k=2)
    print("relative steenrod[5]:", steenrod_rel[5])

    # Absolute convention (opt-in): tuples are (birth, death) with birth <=
    # death (death = +inf for essential), and the index is the absolute image
    # degree.  Sq^2 detects the non-trivial H^4 class; the bar lives in
    # steenrod[4].
    _, steenrod_abs = st.compute_steenrod_barcodes(k=2, absolute=True)
    print("absolute steenrod[4]:", steenrod_abs[4])

.. testoutput::

    relative steenrod[5]: [(0.0, 1.0)]
    absolute steenrod[4]: [(0.0, 1.0)]
