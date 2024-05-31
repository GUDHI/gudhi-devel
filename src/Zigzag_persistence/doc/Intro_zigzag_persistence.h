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
 * We refer to the introduction page \ref persistent_cohomology for persistent (co)homology for an introduction
 * to the topic.
 * Zigzag persistence is a generalization of the latter. While standard persistence only allows to grow the filtered
 * complex by adding simplices, zigzag persistence also allows removals. Hence the name "zigzag", as the module
 * diagram will have arrows alterning between forward and backward.
 *
 * The implementation is based on the algorithm introduced in \cite zigzag.
 *
 * \subsection zigzaginterface Stream-like interface
 * 
 * As removals are possible in zigzag filtration, the maximal size of the complex does not depend on the length of the 
 * filtration anymore. This makes it possible to build very long fine tuned filtrations with relatively small complexes
 * which can be processed without overreaching memory space. For this purpose, it is possible to feed the module with
 * information about the filtration "on the fly" to avoid loading the whole filtration at once. Information about the
 * current barcode can be retrieved between any steps.
 * 
 * \subsection zigzagexamples Examples
 * 
 * Here is a list of zigzag persistence examples :
 * \li \gudhi_example_link{Zigzag_persistence,example_simple_zigzag_filtration.cpp} - A simple example to showcase how
 * to use the \ref Zigzag_persistence class.
 *
 * \li \gudhi_example_link{Zigzag_persistence,example_zzfiltration_from_file.cpp} - An example of a "stream-like" usage
 * by reading of the filtration from a file.
 * 
 * @}
 */
}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // DOC_ZIGZAG_PERSISTENCE_INTRO_ZIGZAG_PERSISTENCE_H_
