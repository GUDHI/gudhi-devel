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
 * \author    Hannah Schreiber, Clément Maria
 *
 * Set of classes allowing addition and multiplication, as well as inverse computation, in \f$ \mathbb{F}_p \f$ fields,
 * with \f$ p \f$ some prime number, or in multi-fields as defined in \cite boissonnat:hal-00922572.
 *
 * There are two types of classes:
 * - those defining directly a field element, allowing to use them as any integer: the operators are overwritten such
 *   that the calculation is done in the field. For example, if \f$ e = 2 \f$ is an instanciation of an
 *   \f$ \mathbb{F}_3 \f$ element class, then `e + 3` returns `2`,
 * - those only defining the operators of a field or multi-field. They represent a collection of methods taking
 *   one or two integers as input and treating them as elements of the field. For example, if \f$ op \f$ is an
 *   instanciation of a \f$ \mathbb{F}_3 \f$ operator class, `op.add(2, 3)` returns `2`.
 *
 * The field operator classes all respect the @ref persistence_matrix::FieldOperators concept.
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
