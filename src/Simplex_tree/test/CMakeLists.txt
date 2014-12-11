cmake_minimum_required(VERSION 2.6)
project(GUDHITestSimplexTree)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} --coverage")
endif()

add_executable ( simplex_tree_unit_test simplex_tree_unit_test.cpp )
target_link_libraries(simplex_tree_unit_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(simplex_tree_unit_test ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_unit_test)

if (LCOV_PATH)
  # Lcov code coverage of unitary test
  add_test(src/Simplex_tree/lcov/coverage.log ${CMAKE_SOURCE_DIR}/scripts/check_code_coverage.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree)
endif()

if (PYTHON_PATH)
  # Cpplint tests on coding style
  add_test(Simplex_tree.h.cpplint ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi/Simplex_tree.h)
  add_test(indexing_tag.h.cpplint ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi/Simplex_tree/indexing_tag.h)
  add_test(Simplex_tree_iterators.h.cpplint ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi/Simplex_tree/Simplex_tree_iterators.h)
  add_test(Simplex_tree_node_explicit_storage.h.cpplint ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi/Simplex_tree/Simplex_tree_node_explicit_storage.h)
  add_test(Simplex_tree_siblings.h.cpplint ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi/Simplex_tree/Simplex_tree_siblings.h)
endif()