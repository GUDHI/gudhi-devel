cmake_minimum_required(VERSION 2.6)
project(Bottleneck_distance_examples)

if (NOT CGAL_VERSION VERSION_LESS 4.8.1)
  add_executable (bottleneck_basic_example bottleneck_basic_example.cpp)
  add_executable (bottleneck_read_file_example bottleneck_read_file_example.cpp)

  add_test(NAME Bottleneck_distance_example_basic COMMAND $<TARGET_FILE:bottleneck_basic_example>)

  add_test(NAME Bottleneck_read_file_example
      COMMAND $<TARGET_FILE:bottleneck_read_file_example>
      "${CMAKE_SOURCE_DIR}/data/persistence_diagram/first.pers" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/second.pers")

  install(TARGETS bottleneck_read_file_example DESTINATION bin)
  install(TARGETS bottleneck_basic_example DESTINATION bin)
  
endif (NOT CGAL_VERSION VERSION_LESS 4.8.1)
