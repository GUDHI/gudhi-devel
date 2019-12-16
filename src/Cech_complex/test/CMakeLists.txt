cmake_minimum_required(VERSION 2.6)
project(Cech_complex_tests)

include(GUDHI_boost_test)

add_executable ( Cech_complex_test_unit test_cech_complex.cpp )
if (TBB_FOUND)
  target_link_libraries(Cech_complex_test_unit ${TBB_LIBRARIES})
endif()

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

gudhi_add_boost_test(Cech_complex_test_unit)
