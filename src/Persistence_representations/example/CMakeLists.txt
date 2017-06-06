cmake_minimum_required(VERSION 2.6)
project(Persistence_representations_example)

file(COPY "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

add_executable ( Persistence_representations_example_landscape_on_grid persistence_landscape_on_grid.cpp )
target_link_libraries(Persistence_representations_example_landscape_on_grid ${Boost_SYSTEM_LIBRARY})
add_test(NAME Persistence_representations_example_landscape_on_grid
    COMMAND $<TARGET_FILE:Persistence_representations_example_landscape_on_grid>)

add_executable ( Persistence_representations_example_landscape persistence_landscape.cpp )
target_link_libraries(Persistence_representations_example_landscape ${Boost_SYSTEM_LIBRARY})
add_test(NAME Persistence_representations_example_landscape
    COMMAND $<TARGET_FILE:Persistence_representations_example_landscape>)

add_executable ( Persistence_representations_example_intervals persistence_intervals.cpp )
target_link_libraries(Persistence_representations_example_intervals ${Boost_SYSTEM_LIBRARY})
add_test(NAME Persistence_representations_example_intervals
    COMMAND $<TARGET_FILE:Persistence_representations_example_intervals>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( Persistence_representations_example_vectors persistence_vectors.cpp )
target_link_libraries(Persistence_representations_example_vectors ${Boost_SYSTEM_LIBRARY})
add_test(NAME Persistence_representations_example_vectors
    COMMAND $<TARGET_FILE:Persistence_representations_example_vectors>)

add_executable ( Persistence_representations_example_heat_maps persistence_heat_maps.cpp )
target_link_libraries(Persistence_representations_example_heat_maps ${Boost_SYSTEM_LIBRARY})
add_test(NAME Persistence_representations_example_heat_maps
    COMMAND $<TARGET_FILE:Persistence_representations_example_heat_maps>)




