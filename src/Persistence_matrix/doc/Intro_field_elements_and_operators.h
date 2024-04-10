/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_PERSISTENCE_MATRIX_INTRO_FIELDS_H_
#define DOC_PERSISTENCE_MATRIX_INTRO_FIELDS_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
namespace persistence_fields {

/** \defgroup persistence_fields Persistence Fields
 * @{
 * \author    Hannah Schreiber
 *
 * The module provides a data structure for matrices, in particular thought for matrices representing filtered complexes
 * and used as backend for persistence algorithms, such at persistent homology, \ref persistent_cohomology,
 * or zigzag. TODO: refs.
 *
 * The structure provides several functionnalities which can be enabled or disabled through a structure following
 * the @ref PersistenceMatrixOptions concept which is then used as template parameter.
 * The main functionnalities are:
 * @li column and row access,
 * @li column addition and scalar multiplication,
 * @li removal of maximal faces while maintaining a valid reduced @ref boundarymatrix "boundary matrix" or compatible chain complex base
 * and a valid barcode with respect to the new filtration,
 * @li computation of persistent homology (but note that if the barcode is your only necessity, using the
 * \ref persistent_cohomology module is often more performant),
 * @li computation of representative cycles for the cycle classes,
 * @li swapping of two consecutive faces in a filtration (cf. vineyards) while maintaining a valid reduced boundary
 * matrix or compatible chain complex base and a valid barcode with respect to the new filtration,
 * 
 * 
 * \subsection fieldsexamples Examples
 * 
 * Here is a list of examples using the module: TODO:
 * \li \gudhi_example_link{Zigzag_persistence,example_simple_zigzag_filtration.cpp} - A simple example to showcase how
 * to use the \ref Zigzag_persistence class.
 *
 * \li \gudhi_example_link{Zigzag_persistence,example_zzfiltration_from_file.cpp} - An example of a "stream-like" usage
 * by reading of the filtration from a file.
 * 
 * @}
 */
}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // DOC_PERSISTENCE_MATRIX_INTRO_FIELDS_H_
