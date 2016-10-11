cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)

add_executable ( persistence_landscape_on_grid persistence_landscape_on_grid.cpp )
target_link_libraries(persistence_landscape_on_grid ${Boost_SYSTEM_LIBRARY})

add_executable ( persistence_landscape persistence_landscape.cpp )
target_link_libraries(persistence_landscape ${Boost_SYSTEM_LIBRARY})

add_executable ( persistence_intervals persistence_intervals.cpp )
target_link_libraries(persistence_intervals ${Boost_SYSTEM_LIBRARY})

add_executable ( vector_representation vector_representation.cpp )
target_link_libraries(vector_representation ${Boost_SYSTEM_LIBRARY})

add_executable ( persistence_heat_maps persistence_heat_maps.cpp )
target_link_libraries(persistence_heat_maps ${Boost_SYSTEM_LIBRARY})




