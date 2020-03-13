project(Simplex_tree_tests)

include(GUDHI_boost_test)

# Do not forget to copy test files in current binary dir
file(COPY "simplex_tree_for_unit_test.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_executable ( Simplex_tree_test_unit simplex_tree_unit_test.cpp )
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_test_unit ${TBB_LIBRARIES})
endif()
gudhi_add_boost_test(Simplex_tree_test_unit)

add_executable ( Simplex_tree_remove_test_unit simplex_tree_remove_unit_test.cpp )
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_remove_test_unit ${TBB_LIBRARIES})
endif()
gudhi_add_boost_test(Simplex_tree_remove_test_unit)

add_executable ( Simplex_tree_iostream_operator_test_unit simplex_tree_iostream_operator_unit_test.cpp )
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_iostream_operator_test_unit ${TBB_LIBRARIES})
endif()
gudhi_add_boost_test(Simplex_tree_iostream_operator_test_unit)

add_executable ( Simplex_tree_ctor_and_move_test_unit simplex_tree_ctor_and_move_unit_test.cpp )
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_ctor_and_move_test_unit ${TBB_LIBRARIES})
endif()
gudhi_add_boost_test(Simplex_tree_ctor_and_move_test_unit)

add_executable ( Simplex_tree_make_filtration_non_decreasing_test_unit simplex_tree_make_filtration_non_decreasing_unit_test.cpp )
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_make_filtration_non_decreasing_test_unit ${TBB_LIBRARIES})
endif()
gudhi_add_boost_test(Simplex_tree_make_filtration_non_decreasing_test_unit)
