cmake_minimum_required(VERSION 2.6)
project(Persistence_representations_example)

file(COPY "simple_diagram.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_executable ( persistence_landscape_on_grid persistence_landscape_on_grid.cpp )
target_link_libraries(persistence_landscape_on_grid ${Boost_SYSTEM_LIBRARY})
add_test ( persistence_landscape_on_grid ${CMAKE_CURRENT_BINARY_DIR}/persistence_landscape_on_grid )

add_executable ( persistence_landscape persistence_landscape.cpp )
target_link_libraries(persistence_landscape ${Boost_SYSTEM_LIBRARY})
add_test ( persistence_landscape ${CMAKE_CURRENT_BINARY_DIR}/persistence_landscape )

add_executable ( persistence_intervals persistence_intervals.cpp )
target_link_libraries(persistence_intervals ${Boost_SYSTEM_LIBRARY})
add_test ( persistence_intervals ${CMAKE_CURRENT_BINARY_DIR}/persistence_intervals simple_diagram.txt )

add_executable ( persistence_vectors persistence_vectors.cpp )
target_link_libraries(persistence_vectors ${Boost_SYSTEM_LIBRARY})
add_test ( persistence_vectors ${CMAKE_CURRENT_BINARY_DIR}/persistence_vectors )

add_executable ( persistence_heat_maps persistence_heat_maps.cpp )
target_link_libraries(persistence_heat_maps ${Boost_SYSTEM_LIBRARY})
add_test ( persistence_heat_maps ${CMAKE_CURRENT_BINARY_DIR}/persistence_heat_maps )




