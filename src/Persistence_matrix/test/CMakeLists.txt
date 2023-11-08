project(Persistence_matrix_tests)

include(GUDHI_boost_test)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_LIBCPP_DEBUG=1")

add_executable ( Persistence_matrix_column_tests_base Persistence_matrix_column_tests_base.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_column_tests_base)

add_executable ( Persistence_matrix_column_tests_boundary Persistence_matrix_column_tests_boundary.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_column_tests_boundary)

add_executable ( Persistence_matrix_column_tests_chain Persistence_matrix_column_tests_chain.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_column_tests_chain)

add_executable ( Persistence_matrix_matrix_tests_z2_base Persistence_matrix_matrix_tests_z2_base.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_z2_base)

add_executable ( Persistence_matrix_matrix_tests_zp_base Persistence_matrix_matrix_tests_zp_base.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_zp_base)

add_executable ( Persistence_matrix_matrix_tests_z2_compression Persistence_matrix_matrix_tests_z2_compression.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_z2_compression)

add_executable ( Persistence_matrix_matrix_tests_zp_compression Persistence_matrix_matrix_tests_zp_compression.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_zp_compression)

add_executable ( Persistence_matrix_matrix_tests_z2_boundary Persistence_matrix_matrix_tests_z2_boundary.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_z2_boundary)

add_executable ( Persistence_matrix_matrix_tests_zp_boundary Persistence_matrix_matrix_tests_zp_boundary.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_zp_boundary)

add_executable ( Persistence_matrix_matrix_tests_z2_ru_vine Persistence_matrix_matrix_tests_z2_ru_vine.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_z2_ru_vine)

add_executable ( Persistence_matrix_matrix_tests_z2_ru_rep Persistence_matrix_matrix_tests_z2_ru_rep.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_z2_ru_rep)

add_executable ( Persistence_matrix_matrix_tests_zp_ru_rep Persistence_matrix_matrix_tests_zp_ru_rep.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_zp_ru_rep)

add_executable ( Persistence_matrix_matrix_tests_z2_chain_vine Persistence_matrix_matrix_tests_z2_chain_vine.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_z2_chain_vine)

add_executable ( Persistence_matrix_matrix_tests_z2_chain_rep Persistence_matrix_matrix_tests_z2_chain_rep.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_z2_chain_rep)

add_executable ( Persistence_matrix_matrix_tests_zp_chain_rep Persistence_matrix_matrix_tests_zp_chain_rep.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_zp_chain_rep)

add_executable ( Persistence_matrix_matrix_tests_z2_chain_barcode Persistence_matrix_matrix_tests_z2_chain_barcode.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_z2_chain_barcode)

add_executable ( Persistence_matrix_matrix_tests_zp_chain_barcode Persistence_matrix_matrix_tests_zp_chain_barcode.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_matrix_tests_zp_chain_barcode)

add_executable ( Persistence_matrix_field_tests Persistence_matrix_field_tests.cpp) #[[${CMAKE_SOURCE_DIR}/debug.cpp)]]
gudhi_add_boost_test(Persistence_matrix_field_tests)

