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

In the **absolute** convention (``absolute=True``, the default) bars satisfy
``birth ≤ death`` and the index is the true cohomological degree of the class.
In the **relative** convention (``absolute=False``) finite bars have ``birth > death``
and the index is one greater than the cohomological degree.

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

    # Absolute convention (default): birth <= death, index = true cohomological degree.
    # Sq^2 detects the non-trivial H^4 class; the bar lives in steenrod[4].
    _, steenrod_abs = st.compute_steenrod_barcodes(k=2)
    print("absolute steenrod[4]:", steenrod_abs[4])

    # Relative convention: birth > death, index shifted by +1.
    # The same bar appears in steenrod[5] with birth and death swapped.
    _, steenrod_rel = st.compute_steenrod_barcodes(k=2, absolute=False)
    print("relative steenrod[5]:", steenrod_rel[5])

.. testoutput::

    absolute steenrod[4]: [(0.0, 1.0)]
    relative steenrod[5]: [(1.0, 0.0)]
