/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_VINEYARD_INTRO_H_
#define DOC_VINEYARD_INTRO_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
namespace vineyard {

/** \defgroup vineyard Vineyard
 * @{
 * \author    Hannah Schreiber
 *
 * \section vineyardintro Vineyard
 *
 * From a filtration of the form
 * \f$ \mathcal{K}_0 \rightarrow \mathcal{K}_1 \rightarrow \cdots \rightarrow \mathcal{K}_n \f$ we can compute the
 * persistence diagram. For example, lets take two adjacent triangles and add they faces one after the other:
 * \code{.cpp}
      Gudhi::Simplex_tree<> st;

      st.insert_simplex({0}, 0);
      st.insert_simplex({1}, 1);
      st.insert_simplex({2}, 2);
      st.insert_simplex({3}, 3);
      st.insert_simplex({0, 3}, 4);
      st.insert_simplex({0, 2}, 5);
      st.insert_simplex({1, 2}, 6);
      st.insert_simplex({2, 3}, 7);
      st.insert_simplex({0, 1}, 8);
      st.insert_simplex({0, 1, 2}, 9);
      st.insert_simplex({0, 2, 3}, 10);
 * \endcode
 * This will result in the following diagram:
 * \image html "vy_intro_ex1.png"
 * But now we could decide to let \f$ \{2, 3\} \f$ appear before \f$ \{1, 2\} \f$ and the resulting filtration would
 * still be valid (any face appears before its cofaces), but the first 1-cycle would be born one step earlier:
 * \image html "vy_intro_ex2.png"
 * Note that the number of bars does not change as the final complex remains the same. And because persistence diagrams
 * are stable under the right distance functions, the points in the diagram will not move a lot. Therefore we can
 * "stack" those diagrams and obtain lines which are tracing the changes of the barcode under those "swap operations".
 * One of those lines is called a **vine** and the set of all the vines of the same filtered complex is called a
 * **vineyard**.
 *
 * \subsection vyimpldetails Implementation details
 *
 * To ensure the right matching of a point in a persistence diagram with the next diagram and to avoid recomputing
 * all persistence pairs from scratch, the first layer of the vineyard is computed like any other diagram, but the
 * remaining layers are computed by updating the columns and rows corresponding to the two swapped cells in the
 * underlying matrix from the initial computation. The implementation is based on @cite vineyards and @cite zigzag.
 *
 * The vineyards are build by the @ref Vineyard_builder class, which uses the @ref Vineyard_base class to compute
 * the updates. The resulting vines can either be represented as flat tensors, or for more convenience, in a @ref Vine
 * class.
 *
 * The update methods can take filtration values which differs completely from the precedent filtration. The sequence
 * of swaps will be automatically deduced from them, but the resulting vineyard will only show the start and end state
 * of the points, not the intermediate swaps.
 *
 * In addition to the vineyard, @ref Vineyard_builder can also retrieve a **representative cycle** for each point in
 * the diagram between each update.
 *
 * Note also the existence of the helper method @ref build_boundary_matrix_from_complex which constructs the right
 * input for @ref Vineyard_builder from a filtered complex such as @ref Gudhi::Simplex_tree and
 * @ref Gudhi::cubical_complex::Bitmap_cubical_complex.
 * 
 * \section vineyardexamples Examples
 *
 * \subsection vyminusage Minimalistic example
 *
 * A simple example to showcase how to use the @ref Vineyard_builder.
 * <details open>
 *    \code{.cpp}
 *      // `true` to save non-trivial representative cycles and `1` to only store 1-cycles
        Gudhi::vineyard::Vineyard_builder<double> vyb(true, 1);

        // underlying complex
        std::vector<std::vector<int>> boundaries = {
            {}, {}, {}, {}, 
            {0, 3}, {0, 2}, {1, 2}, {2, 3}, {0, 1}, 
            {4, 5, 7}, {5, 6, 8}
        };
        std::vector<int> dimensions = {0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2};
        // filtrations values at each step
        std::vector<double> f1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        std::vector<double> f2 = {0, 1, 2, 3, 4, 7, 6, 5, 8, 9, 10};
        std::vector<double> f3 = {0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 10};

        // creates first layer of the vineyard
        vyb.initialize(boundaries, dimensions, f1);
        // if rep 1-cycles are needed
        const auto& cycles_step0 = vyb.get_latest_representative_cycles();
        // do something with `cycles_step0`

        // creates second layer of the vineyard
        vyb.update(f2);
        const auto& cycles_step1 = vyb.get_latest_representative_cycles();
        // do something with `cycles_step1`

        // creates third layer of the vineyard
        vyb.update(f3);
        const auto& cycles_step2 = vyb.get_latest_representative_cycles();
        // do something with `cycles_step2`

        const auto& vineyard = vyb.get_current_vineyard();
        // do something with `vineyard`
        // for example print vines of dimension 1
        for (const auto& vine : vineyard) {
          auto dim = vine.get_dimension();
          if (dim == 1) {
            for (const std::array<double, 2>& pair : vine.get_pairs()) {
              std::cout << "(" << pair[0] << ", " << pair[1] << ") ";
            }
            std::cout << "\n";
          }
        }
        std::cout << "\n";
 *    \endcode
 * </details>
 *
 * \subsection vyexamples More elaborate example
 *
 * Example of how to use the @ref Vineyard_builder to construct the vineyard and representative 1-cycles from point
 * clouds using rips filtrations. For more details, see directly the example file
 * \gudhi_example_link{Vineyard,example_simple_vineyard.cpp}.
 * <details open>
 *   @dontinclude example_simple_vineyard.cpp
 *   @skip int main()
 *   @until }
 * </details>
 * 
 * @}
 */
}  // namespace vineyard
}  // namespace Gudhi

#endif  // DOC_VINEYARD_INTRO_H_
