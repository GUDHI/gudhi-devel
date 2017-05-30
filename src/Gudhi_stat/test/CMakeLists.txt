cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT_CCUT)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
endif()

# copy data directory for tests purpose. QUESTION -- DOEST IT WORK IF WE CREATE A NEW FOLDER FOR EACH TEST? i THINK THAT WE SHOULD MAKE THS COPYING TO EACH TEST FOLDER. 
file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
add_executable ( persistence_intervals_test persistence_intervals_test.cpp )
target_link_libraries(persistence_intervals_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})


# Unitary tests
add_test(NAME persistence_intervals_test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/persistence_intervals_test
  # XML format for Jenkins xUnit plugin 
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/persistence_intervals_UT.xml --log_level=test_suite --report_level=no)

add_executable ( vector_representation_test vector_representation_test.cpp )
target_link_libraries(vector_representation_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME vector_representation_test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/vector_representation_test
  # XML format for Jenkins xUnit plugin 
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/vector_representation_UT.xml --log_level=test_suite --report_level=no)

add_executable ( persistence_lanscapes_test persistence_lanscapes_test.cpp )
target_link_libraries(persistence_lanscapes_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME persistence_lanscapes_test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/persistence_lanscapes_test
  # XML format for Jenkins xUnit plugin 
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/persistence_lanscapes_UT.xml --log_level=test_suite --report_level=no)

add_executable ( persistence_lanscapes_on_grid_test persistence_lanscapes_on_grid_test.cpp )
target_link_libraries(persistence_lanscapes_on_grid_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME persistence_lanscapes_on_grid_test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/persistence_lanscapes_on_grid_test
  # XML format for Jenkins xUnit plugin 
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/persistence_lanscapes_on_grid_UT.xml --log_level=test_suite --report_level=no)

add_executable ( persistence_heat_maps_test persistence_heat_maps_test.cpp )
target_link_libraries(persistence_heat_maps_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME persistence_heat_maps_test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/persistence_heat_maps_test
  # XML format for Jenkins xUnit plugin
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/persistence_heat_maps_UT.xml --log_level=test_suite --report_level=no)

add_executable ( read_persistence_from_file_test read_persistence_from_file_test.cpp )
target_link_libraries(read_persistence_from_file_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME read_persistence_from_file_test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/read_persistence_from_file_test
  # XML format for Jenkins xUnit plugin
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/read_persistence_from_file_UT.xml --log_level=test_suite --report_level=no)



if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
add_executable ( persistence_intervals_with_distances_test persistence_intervals_with_distances_test.cpp )
target_link_libraries(persistence_intervals_with_distances_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
target_link_libraries(persistence_intervals_with_distances_test ${TBB_LIBRARIES})
endif(TBB_FOUND)
# Unitary tests
add_test(NAME persistence_intervals_with_distances_test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/persistence_intervals_with_distances_test
  # XML format for Jenkins xUnit plugin
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/persistence_intervals_with_distances_test_UT.xml --log_level=test_suite --report_level=no)
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)



