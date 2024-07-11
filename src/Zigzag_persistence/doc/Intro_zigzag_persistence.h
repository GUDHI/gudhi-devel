/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_ZIGZAG_PERSISTENCE_INTRO_ZIGZAG_PERSISTENCE_H_
#define DOC_ZIGZAG_PERSISTENCE_INTRO_ZIGZAG_PERSISTENCE_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
namespace zigzag_persistence {

/** \defgroup zigzag_persistence Zigzag Persistence
 * @{
 * \author    Cl&eacute;ment Maria, Hannah Schreiber
 *
 * \section zigzagintro Zigzag Persistence
 *
 * We refer to the introduction page \ref persistent_cohomology for persistent (co)homology for an introduction
 * to the topic.
 * Zigzag persistence is a generalization of the latter. While standard persistence only allows to grow the filtered
 * complex by adding faces, zigzag persistence also allows removals. Hence the name "zigzag", as the module
 * diagram will have arrows alternating between forward and backward.
 *
 * The module is partitioned in two types of classes: filtered by filtration values and filtered by the atomic
 * operations.
 * - There is one atomic class:
 * @ref Zigzag_persistence. It computes the persistence by considering only the index of an atomic operations in the
 * filtration and not its possibly associated filtration value. For example, if a cycle is born at operation number 6
 * and dies at operation number 7, it will output a bar starting at 6 and ending at 7, even if both operations have
 * the same filtration value in the zigzag filtration and therefore the "real" bar has length 0.
 * - There are two filtered classes: @ref Filtered_zigzag_persistence and @ref Filtered_zigzag_persistence_with_storage.
 * They are both based on @ref Zigzag_persistence and manage additionally the filtration values which are ignored by 
 * @ref Zigzag_persistence. They automatically translate the operation numbers into their corresponding filtration
 * values. They also have more flexible inputs (the boundaries do not have to be ordered, nor identified continuously
 * from 0). The two classes diverge on the way they manage the memory: @ref Filtered_zigzag_persistence removes
 * systematically all unnecessary information and outputs a pair as soon it is closed, while
 * @ref Filtered_zigzag_persistence_with_storage will store all informations about filtration values and bars until the
 * end and output the pairs only when asked. Depending on the use and the length of the filtration, one will be more
 * efficient than the other and vice versa.
 *
 * The implementation is based on the algorithm introduced in \cite zigzag.
 *
 * \subsection zigzaginterface Stream-like interface
 * 
 * As removals are possible in zigzag filtration, the maximal size of the complex does not depend on the length of the 
 * filtration anymore. This makes it possible to build very long fine tuned filtrations with relatively small complexes
 * which can be processed without overreaching memory space. For this purpose, it is possible to feed the module with
 * information about the filtration "on the fly" to avoid loading the whole filtration at once. Information about the
 * current barcode can be retrieved between any steps via callback methods.
 * 
 * \subsection zigzagexamples Examples
 * 
 * Here is a list of zigzag persistence examples :
 * \li \gudhi_example_link{Zigzag_persistence,example_simple_zigzag_filtration.cpp} - A simple example to showcase how
 * to use the @ref Filtered_zigzag_persistence_with_storage class.
 *
 * \li \gudhi_example_link{Zigzag_persistence,example_zzfiltration_from_file.cpp} - An example of a "stream-like" usage
 * with @ref Filtered_zigzag_persistence by reading off the filtration from a file.
 * 
 * @}
 */
}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // DOC_ZIGZAG_PERSISTENCE_INTRO_ZIGZAG_PERSISTENCE_H_
