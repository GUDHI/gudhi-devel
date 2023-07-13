project(Zigzag_persistence_tests)

include(GUDHI_boost_test)

add_executable ( Zigzag_persistence_unit_test zigzag_persistence_unit_test.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(Zigzag_persistence_unit_test TBB::tbb)
endif()

# Do not forget to copy test results files in current binary dir
#file(COPY "${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/test/simplex_tree_file_for_unit_test.txt"
#     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# Unitary tests
gudhi_add_boost_test(Zigzag_persistence_unit_test)

