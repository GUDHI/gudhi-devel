#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simplex_tree_remove"
#include <boost/test/unit_test.hpp>

//  ^
// /!\ Nothing else from Simplex_tree shall be included to test includes are well defined.
#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

struct MyOptions : Simplex_tree_options_full_featured {
  // Not doing persistence, so we don't need those
  static const bool store_key = false;
  static const bool store_filtration = false;
  // I have few vertices
  typedef short Vertex_handle;
};

using Mini_stree = Simplex_tree<MyOptions>;
using Stree = Simplex_tree<>;

BOOST_AUTO_TEST_CASE(remove_maximal_simplex) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "REMOVE MAXIMAL SIMPLEX" << std::endl;

  Mini_stree st;

  st.insert_simplex_and_subfaces({0, 1, 6, 7});
  st.insert_simplex_and_subfaces({3, 4, 5});

  // Constructs a copy at this state for further test purpose
  Mini_stree st_pruned = st;

  st.insert_simplex_and_subfaces({3, 0});
  st.insert_simplex_and_subfaces({2, 1, 0});

  // Constructs a copy at this state for further test purpose
  Mini_stree st_complete = st;
  // st_complete and st:
  //    1   6
  //    o---o
  //   /X\7/
  //  o---o---o---o
  //  2   0   3\X/4
  //            o
  //            5
  // st_pruned:
  //    1   6
  //    o---o
  //     \7/
  //      o   o---o
  //      0   3\X/4
  //            o
  //            5

#ifdef GUDHI_DEBUG
  std::cout << "Check exception throw in debug mode" << std::endl;
  // throw excpt because sh has children
  BOOST_CHECK_THROW (st.remove_maximal_simplex(st.find({0, 1, 6})), std::invalid_argument);
  BOOST_CHECK_THROW (st.remove_maximal_simplex(st.find({3})), std::invalid_argument);
  BOOST_CHECK(st == st_complete);
#endif
  std::cout << "st.remove_maximal_simplex({0, 2})" << std::endl;
  st.remove_maximal_simplex(st.find({0, 2}));
  std::cout << "st.remove_maximal_simplex({0, 1, 2})" << std::endl;
  st.remove_maximal_simplex(st.find({0, 1, 2}));
  std::cout << "st.remove_maximal_simplex({1, 2})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 2}));
  std::cout << "st.remove_maximal_simplex({2})" << std::endl;
  st.remove_maximal_simplex(st.find({2}));
  std::cout << "st.remove_maximal_simplex({3})" << std::endl;
  st.remove_maximal_simplex(st.find({0, 3}));
  
  BOOST_CHECK(st == st_pruned);
  // Remove all, but as the simplex tree is not storing filtration, there is no modification
  st.prune_above_filtration(0.0);
  BOOST_CHECK(st == st_pruned);
  
  Mini_stree st_wo_seven;

  st_wo_seven.insert_simplex_and_subfaces({0, 1, 6});
  st_wo_seven.insert_simplex_and_subfaces({3, 4, 5});
  // st_wo_seven:
  //    1   6
  //    o---o
  //     \X/
  //      o   o---o
  //      0   3\X/4
  //            o
  //            5

  // Remove all 7 to test the both remove_maximal_simplex cases (when _members is empty or not)
  std::cout << "st.remove_maximal_simplex({0, 1, 6, 7})" << std::endl;
  st.remove_maximal_simplex(st.find({0, 1, 6, 7}));
  std::cout << "st.remove_maximal_simplex({0, 1, 7})" << std::endl;
  st.remove_maximal_simplex(st.find({0, 1, 7}));
  std::cout << "st.remove_maximal_simplex({0, 6, 7})" << std::endl;
  st.remove_maximal_simplex(st.find({0, 6, 7}));
  std::cout << "st.remove_maximal_simplex({0, 7})" << std::endl;
  st.remove_maximal_simplex(st.find({0, 7}));
  std::cout << "st.remove_maximal_simplex({1, 6, 7})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 6, 7}));
  std::cout << "st.remove_maximal_simplex({1, 7})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 7}));
  std::cout << "st.remove_maximal_simplex({6, 7})" << std::endl;
  st.remove_maximal_simplex(st.find({6, 7}));
  std::cout << "st.remove_maximal_simplex({7})" << std::endl;
  st.remove_maximal_simplex(st.find({7}));
  
  std::cout << "st.dimension()=" << st.dimension() << " | st_wo_seven.dimension()=" << st_wo_seven.dimension() << std::endl;
  BOOST_CHECK(st == st_wo_seven);
}

BOOST_AUTO_TEST_CASE(auto_dimension_set) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "DIMENSION ON REMOVE MAXIMAL SIMPLEX" << std::endl;

  Mini_stree st;

  st.insert_simplex_and_subfaces({0, 1, 2});
  st.insert_simplex_and_subfaces({0, 1, 3});
  st.insert_simplex_and_subfaces({1, 2, 3, 4});
  st.insert_simplex_and_subfaces({1, 2, 3, 5});
  st.insert_simplex_and_subfaces({6, 7, 8, 9});
  st.insert_simplex_and_subfaces({6, 7, 8, 10});

  BOOST_CHECK(st.dimension() == 3);

  std::cout << "st.remove_maximal_simplex({6, 7, 8, 10})" << std::endl;
  st.remove_maximal_simplex(st.find({6, 7, 8, 10}));
  BOOST_CHECK(st.dimension() == 3);

  std::cout << "st.remove_maximal_simplex({6, 7, 8, 9})" << std::endl;
  st.remove_maximal_simplex(st.find({6, 7, 8, 9}));
  BOOST_CHECK(st.dimension() == 3);

  std::cout << "st.remove_maximal_simplex({1, 2, 3, 4})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 2, 3, 4}));
  BOOST_CHECK(st.dimension() == 3);

  std::cout << "st.remove_maximal_simplex({1, 2, 3, 5})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 2, 3, 5}));
  BOOST_CHECK(st.dimension() == 2);

  std::cout << "st.insert_simplex_and_subfaces({1, 2, 3, 5})" << std::endl;
  st.insert_simplex_and_subfaces({1, 2, 3, 5});
  BOOST_CHECK(st.dimension() == 3);

  std::cout << "st.insert_simplex_and_subfaces({1, 2, 3, 4})" << std::endl;
  st.insert_simplex_and_subfaces({1, 2, 3, 4});
  BOOST_CHECK(st.dimension() == 3);

  std::cout << "st.remove_maximal_simplex({1, 2, 3, 5})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 2, 3, 5}));
  BOOST_CHECK(st.dimension() == 3);

  std::cout << "st.remove_maximal_simplex({1, 2, 3, 4})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 2, 3, 4}));
  BOOST_CHECK(st.dimension() == 2);

  std::cout << "st.insert_simplex_and_subfaces({0, 1, 3, 4})" << std::endl;
  st.insert_simplex_and_subfaces({0, 1, 3, 4});
  BOOST_CHECK(st.dimension() == 3);

  std::cout << "st.remove_maximal_simplex({0, 1, 3, 4})" << std::endl;
  st.remove_maximal_simplex(st.find({0, 1, 3, 4}));
  BOOST_CHECK(st.dimension() == 2);

  std::cout << "st.insert_simplex_and_subfaces({1, 2, 3, 5})" << std::endl;
  st.insert_simplex_and_subfaces({1, 2, 3, 5});
  BOOST_CHECK(st.dimension() == 3);

  std::cout << "st.insert_simplex_and_subfaces({1, 2, 3, 4})" << std::endl;
  st.insert_simplex_and_subfaces({1, 2, 3, 4});
  BOOST_CHECK(st.dimension() == 3);

  // Check you can override the dimension
  // This is a limit test case - shall not happen
  st.set_dimension(1);
  BOOST_CHECK(st.dimension() == 1);

  // Here no siblings is erased - automatic dimension is not launched.
  std::cout << "st.remove_maximal_simplex({1, 2, 3, 4})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 2, 3, 4}));
  BOOST_CHECK(st.dimension() == 1);

  // Here sibling is erased - automatic dimension is launched.
  std::cout << "st.remove_maximal_simplex({1, 2, 3, 5})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 2, 3, 5}));
  BOOST_CHECK(st.dimension() == 2);

  std::cout << "st.insert_simplex_and_subfaces({1, 2, 3, 5})" << std::endl;
  st.insert_simplex_and_subfaces({1, 2, 3, 5});
  BOOST_CHECK(st.dimension() == 3);

  std::cout << "st.insert_simplex_and_subfaces({1, 2, 3, 4})" << std::endl;
  st.insert_simplex_and_subfaces({1, 2, 3, 4});
  BOOST_CHECK(st.dimension() == 3);

  // Check you can override the dimension
  // This is a limit test case - shall not happen
  st.set_dimension(6);
  BOOST_CHECK(st.dimension() == 6);

  // Here no siblings is erased - automatic dimension is not launched.
  std::cout << "st.remove_maximal_simplex({1, 2, 3, 5})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 2, 3, 5}));
  BOOST_CHECK(st.dimension() == 6);

  // Here sibling is erased - automatic dimension is launched but dim is always than te one set.
  std::cout << "st.remove_maximal_simplex({1, 2, 3, 4})" << std::endl;
  st.remove_maximal_simplex(st.find({1, 2, 3, 4}));
  BOOST_CHECK(st.dimension() == 6);

  // Reset with the correct value
  st.set_dimension(2);
  BOOST_CHECK(st.dimension() == 2);

  std::cout << "st.insert_simplex_and_subfaces({0, 1, 2, 3, 4, 5, 6})" << std::endl;
  st.insert_simplex_and_subfaces({0, 1, 2, 3, 4, 5, 6});
  BOOST_CHECK(st.dimension() == 6);

  std::cout << "st.remove_maximal_simplex({0, 1, 2, 3, 4, 5, 6})" << std::endl;
  st.remove_maximal_simplex(st.find({0, 1, 2, 3, 4, 5, 6}));
  BOOST_CHECK(st.dimension() == 5);

}

BOOST_AUTO_TEST_CASE(prune_above_filtration) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "PRUNE ABOVE FILTRATION" << std::endl;

  Stree st;

  st.insert_simplex_and_subfaces({0, 1, 6, 7}, 1.0);
  st.insert_simplex_and_subfaces({3, 4, 5}, 2.0);

  // Constructs a copy at this state for further test purpose
  Stree st_pruned = st;
  st_pruned.initialize_filtration();  // reset

  st.insert_simplex_and_subfaces({3, 0}, 3.0);
  st.insert_simplex_and_subfaces({2, 1, 0}, 4.0);

  // Constructs a copy at this state for further test purpose
  Stree st_complete = st;
  // st_complete and st:
  //    1   6
  //    o---o
  //   /X\7/
  //  o---o---o---o
  //  2   0   3\X/4
  //            o
  //            5
  // st_pruned:
  //    1   6
  //    o---o
  //     \7/
  //      o   o---o
  //      0   3\X/4
  //            o
  //            5

  bool simplex_is_changed = false;
  // Check the no action cases
  // greater than initial filtration value
  simplex_is_changed = st.prune_above_filtration(10.0);
  if (simplex_is_changed)
    st.initialize_filtration();
  BOOST_CHECK(st == st_complete);
  BOOST_CHECK(!simplex_is_changed);
  // equal to initial filtration value
  simplex_is_changed = st.prune_above_filtration(6.0);
  if (simplex_is_changed)
    st.initialize_filtration();
  BOOST_CHECK(st == st_complete);
  BOOST_CHECK(!simplex_is_changed);
  // lower than initial filtration value, but still greater than the maximum filtration value
  simplex_is_changed = st.prune_above_filtration(5.0);
  if (simplex_is_changed)
    st.initialize_filtration();
  BOOST_CHECK(st == st_complete);
  BOOST_CHECK(!simplex_is_changed);

  // Display the Simplex_tree
  std::cout << "The complex contains " << st.num_simplices() << " simplices";
  std::cout << " - dimension " << st.dimension() << std::endl;
  std::cout << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::cout << (int) vertex << " ";
    }
    std::cout << std::endl;
  }

  // Check the pruned cases
  simplex_is_changed = st.prune_above_filtration(2.5);
  if (simplex_is_changed)
    st.initialize_filtration();
  BOOST_CHECK(st == st_pruned);
  BOOST_CHECK(simplex_is_changed);

  // Display the Simplex_tree
  std::cout << "The complex pruned at 2.5 contains " << st.num_simplices() << " simplices";
  std::cout << " - dimension " << st.dimension() << std::endl;

  simplex_is_changed = st.prune_above_filtration(2.0);
  if (simplex_is_changed)
    st.initialize_filtration();
  
  std::cout << "The complex pruned at 2.0 contains " << st.num_simplices() << " simplices";
  std::cout << " - dimension " << st.dimension() << std::endl;

  BOOST_CHECK(st == st_pruned);
  BOOST_CHECK(!simplex_is_changed);

  Stree st_empty;
  simplex_is_changed = st.prune_above_filtration(0.0);
  if (simplex_is_changed)
    st.initialize_filtration();

  // Display the Simplex_tree
  std::cout << "The complex pruned at 0.0 contains " << st.num_simplices() << " simplices";
  std::cout << " - dimension " << st.dimension() << std::endl;

  BOOST_CHECK(st == st_empty);
  BOOST_CHECK(simplex_is_changed);

  // Test case to the limit
  simplex_is_changed = st.prune_above_filtration(-1.0);
  if (simplex_is_changed)
    st.initialize_filtration();
  BOOST_CHECK(st == st_empty);
  BOOST_CHECK(!simplex_is_changed);
}

BOOST_AUTO_TEST_CASE(mini_prune_above_filtration) {
  std::cout << "********************************************************************" << std::endl;
  std::cout << "MINI PRUNE ABOVE FILTRATION" << std::endl;

  Mini_stree st;

  st.insert_simplex_and_subfaces({0, 1, 6, 7});
  st.insert_simplex_and_subfaces({3, 4, 5});
  st.insert_simplex_and_subfaces({3, 0});
  st.insert_simplex_and_subfaces({2, 1, 0});

  // st:
  //    1   6
  //    o---o
  //   /X\7/
  //  o---o---o---o
  //  2   0   3\X/4
  //            o
  //            5

  st.initialize_filtration();
  
  // Display the Simplex_tree
  std::cout << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
  BOOST_CHECK(st.num_simplices() == 27);

  // Test case to the limit - With these options, there is no filtration, which means filtration is 0
  bool simplex_is_changed = st.prune_above_filtration(1.0);
  if (simplex_is_changed)
    st.initialize_filtration();
  // Display the Simplex_tree
  std::cout << "The complex pruned at 1.0 contains " << st.num_simplices() << " simplices" << std::endl;
  BOOST_CHECK(!simplex_is_changed);
  BOOST_CHECK(st.num_simplices() == 27);

  simplex_is_changed = st.prune_above_filtration(0.0);
  if (simplex_is_changed)
    st.initialize_filtration();
  // Display the Simplex_tree
  std::cout << "The complex pruned at 0.0 contains " << st.num_simplices() << " simplices" << std::endl;
  BOOST_CHECK(!simplex_is_changed);
  BOOST_CHECK(st.num_simplices() == 27);

  // Test case to the limit
  simplex_is_changed = st.prune_above_filtration(-1.0);
  if (simplex_is_changed)
    st.initialize_filtration();
  // Display the Simplex_tree
  std::cout << "The complex pruned at -1.0 contains " << st.num_simplices() << " simplices" << std::endl;
  BOOST_CHECK(simplex_is_changed);
  BOOST_CHECK(st.num_simplices() == 0);

  // Display the Simplex_tree
  std::cout << "The complex contains " << st.num_simplices() << " simplices" << std::endl;

}
