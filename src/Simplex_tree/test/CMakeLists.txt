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
add_executable ( TestSimplexTree UnitTestSimplexTree.cpp )
target_link_libraries(TestSimplexTree ${Boost_LIBRARIES})

# Unitary tests
add_test(TestSimplexTree ${CMAKE_CURRENT_BINARY_DIR}/TestSimplexTree)

if (LCOV_PATH)
  # Lcov code coverage of unitary test
  add_test(TestSimplexTreeCov ${CMAKE_SOURCE_DIR}/scripts/check_code_coverage.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree)
endif()

if (PYTHON_PATH)
  # Cpplint tests on coding style
  add_test(CpplintSimplexTree ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi/Simplex_tree.h)
  add_test(CpplintSimplexTreeIndexingTag ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi/Simplex_tree/indexing_tag.h)
  add_test(CpplintSimplexTreeIterator ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi/Simplex_tree/Simplex_tree_iterators.h)
  add_test(CpplintSimplexTreeNode ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi/Simplex_tree/Simplex_tree_node_explicit_storage.h)
  add_test(CpplintSimplexTreeSiblings ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi/Simplex_tree/Simplex_tree_siblings.h)
endif()