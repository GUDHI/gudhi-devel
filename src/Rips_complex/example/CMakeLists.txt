cmake_minimum_required(VERSION 2.6)
project(Rips_complex_examples)

add_executable ( ripsoffreader example_rips_complex_from_off_file.cpp )
target_link_libraries(ripsoffreader  ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(ripsoffreader ${TBB_LIBRARIES})
endif()
