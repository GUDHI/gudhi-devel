/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Anibal M. Medina-Mardones
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_STEENROD_PERSISTENCE_INTRO_STEENROD_PERSISTENCE_H_
#define DOC_STEENROD_PERSISTENCE_INTRO_STEENROD_PERSISTENCE_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace steenrod_persistence {

/** \defgroup steenrod_persistence Persistence Steenrod modules

  \author    Anibal M. Medina-Mardones

  Given a filtered simplicial complex this module computes, over
  \f$\mathbb{F}_2\f$, the ordinary persistence barcode together with the
  Sq\f$^k\f$ Steenrod barcode. The Steenrod barcode is an algebraic-topology
  invariant generalizing the ordinary mod 2 persistence barcode. This
  module is based on the project Steenroder by Lupo, Medina-Mardones,
  and Tauzin.

  \section implementation Implementation overview

  Coefficients are fixed to \f$\mathbb{F}_2\f$ — Steenrod squares are defined
  only over the prime field of characteristic two — so the module exposes no
  field parameter. The pipeline runs in four stages:

    1. Matrix reduction in relative cohomology with the clearing optimisation
       produces an \f$R = DV\f$ decomposition per dimension.
    2. Ordinary persistence pairs and cocycle representatives are read off
       from \f$R\f$ and \f$V\f$.
    3. Sq\f$^k\f$ is applied to every representative using the explicit cup-i
       formulas of \cite medina2023cupifast, producing a Steenrod matrix.
    4. An augmented column reduction of the Steenrod matrix against the
       previous-dimension \f$R\f$ produces the Sq\f$^k\f$ Steenrod barcode.

  The algorithm follows \cite lupo2022persistencesteenrod. Stages 3 and 4 are
  parallelised with OpenMP when available.

  \section userinterface User interface

  <b>Python.</b> The module is exposed on the Python simplex tree as
  ``gudhi.SimplexTree.compute_steenrod_barcodes(k)``. Following the same
  ergonomic as ``compute_persistence``, it returns an object with both the
  ordinary barcode and the Sq\f$^k\f$ Steenrod barcode, indexed by cohomological
  dimension. No prior call to ``compute_persistence`` is required.

  <b>C++.</b> Pass a \ref Filtration_by_dim directly to \ref barcodes to obtain
  a \ref Barcodes_result. There is no C++ method on \c Simplex_tree — this
  follows the same convention as \c Persistent_cohomology, which is also
  invoked through a separate algorithm class rather than a simplex tree method.

  \section references References

    - U. Lupo, A. M. Medina-Mardones, G. Tauzin,
      *Persistence Steenrod modules*,
      J. Appl. Comput. Topol. (2022).
      https://doi.org/10.1007/s41468-022-00093-7
    - A. M. Medina-Mardones,
      *New formulas for cup-i products and fast computation of Steenrod squares*,
      Comput. Geom. (2023).
      https://doi.org/10.1016/j.comgeo.2023.101999

*/
// end of group steenrod_persistence

}  // namespace steenrod_persistence

}  // namespace Gudhi

#endif  // DOC_STEENROD_PERSISTENCE_INTRO_STEENROD_PERSISTENCE_H_
