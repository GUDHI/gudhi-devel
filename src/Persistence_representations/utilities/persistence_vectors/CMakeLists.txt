cmake_minimum_required(VERSION 2.6)
project(Persistence_vectors_utilities)

file(COPY "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

add_executable ( create_persistence_vectors create_persistence_vectors.cpp )
target_link_libraries(create_persistence_vectors ${Boost_SYSTEM_LIBRARY})

# Will create simple_diagram.txt.vect
add_test(NAME create_persistence_vectors COMMAND $<TARGET_FILE:create_persistence_vectors>
    "-1" "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( average_persistence_vectors average_persistence_vectors.cpp )
target_link_libraries(average_persistence_vectors ${Boost_SYSTEM_LIBRARY})

add_test(NAME average_persistence_vectors COMMAND $<TARGET_FILE:average_persistence_vectors>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.vect"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.vect")

add_executable ( compute_distance_of_persistence_vectors compute_distance_of_persistence_vectors.cpp )
target_link_libraries(compute_distance_of_persistence_vectors ${Boost_SYSTEM_LIBRARY})

add_test(NAME compute_distance_of_persistence_vectors COMMAND $<TARGET_FILE:compute_distance_of_persistence_vectors>
    "-1" "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.vect"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.vect")

add_executable ( compute_scalar_product_of_persistence_vectors compute_scalar_product_of_persistence_vectors.cpp )
target_link_libraries(compute_scalar_product_of_persistence_vectors ${Boost_SYSTEM_LIBRARY})

add_test(NAME compute_scalar_product_of_persistence_vectors COMMAND $<TARGET_FILE:compute_scalar_product_of_persistence_vectors>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.vect"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.vect")

add_executable ( plot_persistence_vectors plot_persistence_vectors.cpp )
target_link_libraries(plot_persistence_vectors ${Boost_SYSTEM_LIBRARY})

add_test(NAME plot_persistence_vectors COMMAND $<TARGET_FILE:plot_persistence_vectors>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.vect")
