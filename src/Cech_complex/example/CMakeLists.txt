cmake_minimum_required(VERSION 2.6)
project(Cech_complex_examples)

add_executable ( Cech_complex_example_step_by_step cech_complex_step_by_step.cpp )
target_link_libraries(Cech_complex_example_step_by_step ${Boost_PROGRAM_OPTIONS_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Cech_complex_example_step_by_step ${TBB_LIBRARIES})
endif()

add_executable ( Cech_complex_example_from_points cech_complex_example_from_points.cpp)
if (TBB_FOUND)
  target_link_libraries(Cech_complex_example_from_points ${TBB_LIBRARIES})
endif()
add_test(NAME Cech_complex_example_from_points COMMAND $<TARGET_FILE:Cech_complex_example_from_points>)
