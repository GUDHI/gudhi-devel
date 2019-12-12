project(Skeleton_blocker_tests)

include(GUDHI_boost_test)

add_executable ( Skeleton_blocker_test_unit test_skeleton_blocker_complex.cpp )
add_executable ( Skeleton_blocker_test_geometric_complex test_skeleton_blocker_geometric_complex.cpp )
add_executable ( Skeleton_blocker_test_simplifiable test_skeleton_blocker_simplifiable.cpp )

# Do not forget to copy test files in current binary dir
file(COPY "test2.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

gudhi_add_boost_test(Skeleton_blocker_test_unit)
gudhi_add_boost_test(Skeleton_blocker_test_geometric_complex)
gudhi_add_boost_test(Skeleton_blocker_test_simplifiable)
