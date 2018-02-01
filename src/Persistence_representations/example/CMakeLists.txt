cmake_minimum_required(VERSION 2.6)
project(Persistence_representations_example)

add_executable ( Persistence_representations_example_landscape_on_grid persistence_landscape_on_grid.cpp )
add_test(NAME Persistence_representations_example_landscape_on_grid
    COMMAND $<TARGET_FILE:Persistence_representations_example_landscape_on_grid>)

add_executable ( Persistence_representations_example_landscape persistence_landscape.cpp )
add_test(NAME Persistence_representations_example_landscape
    COMMAND $<TARGET_FILE:Persistence_representations_example_landscape>)

add_executable ( Persistence_representations_example_intervals persistence_intervals.cpp )
add_test(NAME Persistence_representations_example_intervals
    COMMAND $<TARGET_FILE:Persistence_representations_example_intervals>
    "${CMAKE_SOURCE_DIR}/data/persistence_diagram/first.pers")

add_executable ( Persistence_representations_example_vectors persistence_vectors.cpp )
add_test(NAME Persistence_representations_example_vectors
    COMMAND $<TARGET_FILE:Persistence_representations_example_vectors>)

add_executable ( Persistence_representations_example_heat_maps persistence_heat_maps.cpp )
add_test(NAME Persistence_representations_example_heat_maps
    COMMAND $<TARGET_FILE:Persistence_representations_example_heat_maps>)




