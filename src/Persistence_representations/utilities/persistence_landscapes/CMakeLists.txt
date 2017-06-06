cmake_minimum_required(VERSION 2.6)
project(Persistence_representations_lanscapes_utilities)

file(COPY "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

add_executable ( create_landscapes create_landscapes.cpp )
target_link_libraries(create_landscapes ${Boost_SYSTEM_LIBRARY})

# Will create simple_diagram.txt.land
add_test(NAME create_landscapes COMMAND $<TARGET_FILE:create_landscapes>
    "-1" "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( average_landscapes average_landscapes.cpp )
target_link_libraries(average_landscapes ${Boost_SYSTEM_LIBRARY})

add_test(NAME average_landscapes COMMAND $<TARGET_FILE:average_landscapes>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( plot_landscapes plot_landscapes.cpp )
target_link_libraries(plot_landscapes ${Boost_SYSTEM_LIBRARY})

add_test(NAME plot_landscapes COMMAND $<TARGET_FILE:plot_landscapes>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( compute_distance_of_landscapes compute_distance_of_landscapes.cpp )
target_link_libraries(compute_distance_of_landscapes ${Boost_SYSTEM_LIBRARY})

add_test(NAME compute_distance_of_landscapes COMMAND $<TARGET_FILE:compute_distance_of_landscapes>
    "1" "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.land"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.land")

add_executable ( compute_scalar_product_of_landscapes compute_scalar_product_of_landscapes.cpp )
target_link_libraries(compute_scalar_product_of_landscapes ${Boost_SYSTEM_LIBRARY})

add_test(NAME compute_scalar_product_of_landscapes COMMAND $<TARGET_FILE:compute_scalar_product_of_landscapes>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.land"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.land")
