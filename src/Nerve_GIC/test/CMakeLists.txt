include(GUDHI_boost_test)

add_executable_with_targets(Nerve_GIC_test_unit test_GIC.cpp TBB::tbb)

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

gudhi_add_boost_test(Nerve_GIC_test_unit)
