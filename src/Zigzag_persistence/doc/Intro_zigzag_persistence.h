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
 * The module consists of the @ref Zigzag_persistence class and two wrappers @ref Filtered_zigzag_persistence and
 * @ref Filtered_zigzag_persistence_with_storage "":
 * - @ref Zigzag_persistence computes the persistence of a sequence of insertions and removals. A face can be inserted
 * or removed one at a time and the returned persistence pairs / bars are indexed on the operation numbers.
 * For example, if a cycle is born at operation number 6 and dies at operation number 7, it will output a bar starting
 * at 6 and ending at 7.
 * - @ref Filtered_zigzag_persistence and @ref Filtered_zigzag_persistence_with_storage are adding the notion of
 * "filtration value" to @ref Zigzag_persistence. At each call, an operation can be associated to a filtration value,
 * which will be used to index the returned bars instead (bars with new length 0 are then ignored). The two classes
 * also have more flexible inputs (the boundaries do not have to be ordered, nor identified continuously
 * from 0). The difference between both classes is on the way they manage the memory: @ref Filtered_zigzag_persistence
 * removes systematically all unnecessary information and outputs a pair as soon it is closed, while
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
 * \section zigzagexamples Examples
 *
 * \subsection zzminusage Minimalistic example of usage
 *
 * ### Includes
 *
 * #### Zigzag_persistence
 * ```
 * #include <gudhi/zigzag_persistence.h>
 * ```
 * #### Filtered_zigzag_persistence and Filtered_zigzag_persistence_with_storage
 * ```
 * #include <gudhi/filtered_zigzag_persistence.h>
 * ```
 *
 * ### Useful aliases
 *
 * ```
 * using Zigzag_persistence = Gudhi::zigzag_persistence::Zigzag_persistence<>;
 * using Filtered_zigzag_persistence = Gudhi::zigzag_persistence::Filtered_zigzag_persistence<>;
 * using Filtered_zigzag_persistence_with_storage = Gudhi::zigzag_persistence::Filtered_zigzag_persistence_with_storage<>;
 *
 * using dimension_type = Zigzag_persistence::dimension_type;
 * using index_type = Zigzag_persistence::index;
 * using filtration_value_type = Filtered_zigzag_persistence::filtration_value;
 * ```
 *
 * ### Construction with default values
 *
 * #### Zigzag_persistence
 * ```
 * //Zigzag_persistence(callback) with for example callback method as a anonymous lambda
 * Zigzag_persistence zp([](dimension_type dim, index_type birth, index_type death) {
 *   std::cout << "[" << dim << "] " << birth << " - " << death << std::endl;
 * });
 * ```
 *
 * #### Filtered_zigzag_persistence
 * ```
 * //Filtered_zigzag_persistence(callback) with for example callback method as a anonymous lambda
 * Filtered_zigzag_persistence zp([](dimension_type dim, filtration_value_type birth, filtration_value_type death) {
 *   std::cout << "[" << dim << "] " << birth << " - " << death << std::endl;
 * });
 * ```
 *
 * #### Filtered_zigzag_persistence_with_storage
 * ```
 * Filtered_zigzag_persistence_with_storage zp;
 * ```
 *
 * ### Input of the zigzag sequence/filtration
 *
 * In all cases, it is important that the operations of insertions and removals are made **in the same order**
 * as in the zigzag filtration ones wants to compute the barcode from.
 *
 * #### Zigzag_persistence
 *
 * A face has to be identified in the boundaries by the operation number the face was inserted with in the sequence.
 *
 * ```
 * //inserts vertex 0 -> birth at 0 of 0-cycle
 * zp.insert_face({}, 0);
 * //inserts vertex 1 -> birth at 1 of 0-cycle
 * zp.insert_face({}, 0);
 * //inserts edge 2 = (0,1) -> death at 2 -> outputs (0, 1, 2)
 * zp.insert_face({0, 1}, 1);
 * //inserts vertex 3 -> birth at 3 of 0-cycle
 * zp.insert_face({}, 0);
 * //inserts edge 4 = (0,3) -> death at 4 -> outputs (0, 3, 4)
 * zp.insert_face({0, 3}, 1);
 * //inserts edge 5 = (1,3) -> birth at 5 of 1-cycle
 * zp.insert_face({1, 3}, 1);
 * //removes edge 4 -> death at 6 -> outputs (1, 5, 6)
 * zp.remove_face(4, 1);
 * //removes edge 2 -> birth at 7 of 0-cycle
 * zp.remove_face(2, 1);
 * ```
 *
 * #### Filtered_zigzag_persistence and Filtered_zigzag_persistence_with_storage
 *
 * A face can be identified in the boundaries by any given numerical label, it is just important that the given
 * filtration values are monotonous (ie., either only increasing or only decreasing).
 *
 * ```
 * //inserts vertex 2 at filtration value 0.1 -> birth at 0.1 of 0-cycle
 * zp.insert_face(2, {}, 0, 0.1);
 * //inserts vertex 4 at filtration value 0.1 -> birth at 0.1 of 0-cycle
 * zp.insert_face(4, {}, 0, 0.1);
 * //inserts edge 5 = (2,4) at filtration value 0.3 -> death at 0.3 -> outputs/stores (0, 0.1, 0.3)
 * zp.insert_face(5, {2, 4}, 1, 0.3);
 * //inserts vertex 3 at filtration value 0.4 -> birth at 0.4 of 0-cycle
 * zp.insert_face(3, {}, 0, 0.4);
 * //inserts edge 6 = (2,3) at filtration value 0.4 -> death at 0.4 of the cycle born at 0.4 -> outputs/stores nothing
 * zp.insert_face(6, {2, 3}, 1, 0.4);
 * //inserts edge 9 = (3,4) at filtration value 1.2 -> birth at 1.2 of 1-cycle
 * zp.insert_face(9, {4, 3}, 1, 1.2);
 * //removes edge 6 at filtration value 1.5 -> death at 1.5 -> outputs/stores (1, 1.2, 1.5)
 * zp.remove_face(6, 1, 1.5);
 * //removes edge 5 at filtration value 2.0 -> birth at 2.0 of 0-cycle
 * zp.remove_face(5, 1, 2.0);
 * ```
 *
 * ### Finalizations
 *
 * For Zigzag_persistence and Filtered_zigzag_persistence, only the closed bars where output so far, so
 * the open/infinite bars still need to be retrieved. For Filtered_zigzag_persistence_with_storage, the bars are
 * stored within the class and where not output at all for now.
 *
 * #### Zigzag_persistence
 * ```
 * //retrieve infinite bars, that is cycles which were born but did not die.
 * //in this example, outputs (0, 0) and (0, 7)
 * zp.get_current_infinite_intervals([](dimension_type dim, index_type birth){
 *   std::cout << "[" << dim << "] " << birth << " - inf" << std::endl;
 * });
 * ```
 *
 * #### Filtered_zigzag_persistence
 * ```
 * //retrieve infinite bars, that is cycles which were born but did not die.
 * //in this example, outputs (0, 0.1) and (0, 2.0)
 * zp.get_current_infinite_intervals([](dimension_type dim, filtration_value_type birth){
 *   std::cout << "[" << dim << "] " << birth << " - inf" << std::endl;
 * });
 * ```
 *
 * #### Filtered_zigzag_persistence_with_storage
 * ```
 * //get all bars in a vector
 * auto barcode = zp.get_persistence_diagram();
 *
 * //do something with the vector, e.g., stream out content:
 * for (auto& bar : barcode) {
 *   std::cout << bar << std::endl;
 * }
 * ```
 *
 * \subsection zzexamples More elaborate examples
 * 
 * \li \gudhi_example_link{Zigzag_persistence,example_simple_zigzag_filtration.cpp} - A simple example to showcase how
 * to use the @ref Filtered_zigzag_persistence_with_storage class within an input loop.
 *
 * \li \gudhi_example_link{Zigzag_persistence,example_zzfiltration_from_file.cpp} - An example of a "stream-like" usage
 * with @ref Filtered_zigzag_persistence by reading off the filtration from a file.
 * 
 * @}
 */
}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // DOC_ZIGZAG_PERSISTENCE_INTRO_ZIGZAG_PERSISTENCE_H_
