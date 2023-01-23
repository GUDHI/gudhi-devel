project(Persistence_matrix_tests)

include(GUDHI_boost_test)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_LIBCPP_DEBUG=1")

add_executable ( Persistence_matrix_common_test_unit Persistence_matrix_common_unit_test.cpp )#${CMAKE_SOURCE_DIR}/debug.cpp)
gudhi_add_boost_test(Persistence_matrix_common_test_unit)

add_executable ( Persistence_matrix_specialized_test_unit Persistence_matrix_specialized_unit_test.cpp )#${CMAKE_SOURCE_DIR}/debug.cpp)
gudhi_add_boost_test(Persistence_matrix_specialized_test_unit)

add_executable ( Persistence_matrix_field_test_unit Persistence_matrix_field_unit_test.cpp )#${CMAKE_SOURCE_DIR}/debug.cpp)
gudhi_add_boost_test(Persistence_matrix_field_test_unit)

add_executable ( Persistence_matrix_columns_test Persistence_matrix_columns_test.cpp )#${CMAKE_SOURCE_DIR}/debug.cpp)
gudhi_add_boost_test(Persistence_matrix_columns_test)

add_executable ( Persistence_matrix_compression_test Persistence_matrix_column_compression_unit_test.cpp )#${CMAKE_SOURCE_DIR}/debug.cpp)
gudhi_add_boost_test(Persistence_matrix_compression_test)
