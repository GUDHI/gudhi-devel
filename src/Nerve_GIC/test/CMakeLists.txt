project(Graph_induced_complex_tests)

include(GUDHI_boost_test)

add_executable ( Nerve_GIC_test_unit test_GIC.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(Nerve_GIC_test_unit TBB::tbb)
endif()

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

gudhi_add_boost_test(Nerve_GIC_test_unit)
