cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)

#persitence vectors
add_executable ( create_persistence_vectors create_persistence_vectors.cpp )
target_link_libraries(create_persistence_vectors ${Boost_SYSTEM_LIBRARY})
add_test ( create_persistence_vectors ${CMAKE_CURRENT_BINARY_DIR}/create_persistence_vectors -1 "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt" )

add_executable ( average_persistence_vectors average_persistence_vectors.cpp )
target_link_libraries(average_persistence_vectors ${Boost_SYSTEM_LIBRARY})
add_test ( average_persistence_vectors ${CMAKE_CURRENT_BINARY_DIR}/average_persistence_vectors "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.vect" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.vect" )

add_executable ( compute_distance_of_persistence_vectors compute_distance_of_persistence_vectors.cpp )
target_link_libraries(compute_distance_of_persistence_vectors ${Boost_SYSTEM_LIBRARY})
add_test ( compute_distance_of_persistence_vectors ${CMAKE_CURRENT_BINARY_DIR}/compute_distance_of_persistence_vectors -1 "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.vect" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.vect" )

add_executable ( compute_scalar_product_of_persistence_vectors compute_scalar_product_of_persistence_vectors.cpp )
target_link_libraries(compute_scalar_product_of_persistence_vectors ${Boost_SYSTEM_LIBRARY})
add_test ( compute_scalar_product_of_persistence_vectors ${CMAKE_CURRENT_BINARY_DIR}/compute_scalar_product_of_persistence_vectors "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.vect" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.vect" )

add_executable ( plot_persistence_vectors plot_persistence_vectors.cpp )
target_link_libraries(plot_persistence_vectors ${Boost_SYSTEM_LIBRARY})
add_test ( plot_persistence_vectors ${CMAKE_CURRENT_BINARY_DIR}/plot_persistence_vectors "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.vect" )


