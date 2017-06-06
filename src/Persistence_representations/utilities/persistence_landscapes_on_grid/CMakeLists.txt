cmake_minimum_required(VERSION 2.6)
project(Persistence_representations_lanscapes_on_grid_utilities)

file(COPY "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

add_executable ( create_landscapes_on_grid create_landscapes_on_grid.cpp )
target_link_libraries(create_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})

# Will create simple_diagram.txt.land
add_test(NAME create_landscapes_on_grid COMMAND $<TARGET_FILE:create_landscapes_on_grid>
    "100" "-1" "-1" "-1" "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( average_landscapes_on_grid average_landscapes_on_grid.cpp )
target_link_libraries(average_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})

add_test(NAME average_landscapes_on_grid COMMAND $<TARGET_FILE:average_landscapes_on_grid>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.g_land"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.g_land")

add_executable ( plot_landscapes_on_grid plot_landscapes_on_grid.cpp )
target_link_libraries(plot_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})

add_test(NAME plot_landscapes_on_grid COMMAND $<TARGET_FILE:plot_landscapes_on_grid>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.g_land")

add_executable ( compute_distance_of_landscapes_on_grid compute_distance_of_landscapes_on_grid.cpp )
target_link_libraries(compute_distance_of_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})

add_test(NAME compute_distance_of_landscapes_on_grid COMMAND $<TARGET_FILE:compute_distance_of_landscapes_on_grid>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.g_land"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.g_land")

add_executable ( compute_scalar_product_of_landscapes_on_grid compute_scalar_product_of_landscapes_on_grid.cpp )
target_link_libraries(compute_scalar_product_of_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})

add_test(NAME compute_scalar_product_of_landscapes_on_grid COMMAND $<TARGET_FILE:compute_scalar_product_of_landscapes_on_grid>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.g_land"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt.g_land")
