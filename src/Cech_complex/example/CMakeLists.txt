cmake_minimum_required(VERSION 2.6)
project(Cech_complex_examples)

add_executable ( Cech_complex_example_step_by_step cech_complex_step_by_step.cpp )
target_link_libraries(Cech_complex_example_step_by_step Boost::program_options)
if (TBB_FOUND)
  target_link_libraries(Cech_complex_example_step_by_step ${TBB_LIBRARIES})
endif()
add_test(NAME Cech_complex_utility_from_rips_on_tore_3D COMMAND $<TARGET_FILE:Cech_complex_example_step_by_step>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-r" "0.25" "-d" "3")

add_executable ( Cech_complex_example_from_points cech_complex_example_from_points.cpp)
if (TBB_FOUND)
  target_link_libraries(Cech_complex_example_from_points ${TBB_LIBRARIES})
endif()
add_test(NAME Cech_complex_example_from_points COMMAND $<TARGET_FILE:Cech_complex_example_from_points>)
