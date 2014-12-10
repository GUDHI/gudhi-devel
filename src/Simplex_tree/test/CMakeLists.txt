cmake_minimum_required(VERSION 2.6)
project(GUDHITestSimplexTree)

# NEED TO FIND BOOST NEEDED COMPONENTS TO LINK WITH
find_package(Boost 1.45.0 COMPONENTS system unit_test_framework)
message("Boost_lib = ${Boost_LIBRARIES}")

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} --coverage")
endif()

include_directories(${Boost_INCLUDE_DIRS})
add_executable ( simplex_tree_unit_test simplex_tree_unit_test.cpp )
target_link_libraries(simplex_tree_unit_test ${Boost_LIBRARIES})

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