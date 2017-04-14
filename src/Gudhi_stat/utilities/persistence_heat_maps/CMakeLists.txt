cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)

#persitence heat maps

add_executable ( create_persistence_heat_maps create_persistence_heat_maps.cpp )
target_link_libraries(create_persistence_heat_maps ${Boost_SYSTEM_LIBRARY})
add_test ( create_persistence_heat_maps ${CMAKE_CURRENT_BINARY_DIR}/create_persistence_heat_maps   10 -1 -1 4 -1  ../../example/simple_diagram.txt )

add_executable ( create_pssk create_pssk.cpp )
target_link_libraries(create_pssk ${Boost_SYSTEM_LIBRARY})
add_test ( create_pssk ${CMAKE_CURRENT_BINARY_DIR}/create_pssk  10 -1 -1 4 -1  ../../example/simple_diagram.txt )

add_executable ( create_persistence_heat_maps_weighted_by_distance_from_diagonal create_persistence_heat_maps_weighted_by_distance_from_diagonal.cpp  )
target_link_libraries(create_persistence_heat_maps_weighted_by_distance_from_diagonal ${Boost_SYSTEM_LIBRARY})
add_test ( create_persistence_heat_maps_weighted_by_distance_from_diagonal ${CMAKE_CURRENT_BINARY_DIR}/create_persistence_heat_maps_weighted_by_distance_from_diagonal 10 -1 -1 4 -1 ../../example/simple_diagram.txt )

add_executable ( create_persistence_heat_maps_weighted_by_squared_distance_from_diagonal create_persistence_heat_maps_weighted_by_squared_distance_from_diagonal.cpp )
target_link_libraries(create_persistence_heat_maps_weighted_by_squared_distance_from_diagonal ${Boost_SYSTEM_LIBRARY})
add_test ( create_persistence_heat_maps_weighted_by_squared_distance_from_diagonal ${CMAKE_CURRENT_BINARY_DIR}/create_persistence_heat_maps_weighted_by_squared_distance_from_diagonal 10 -1 -1 4 -1 ../../example/simple_diagram.txt )

add_executable ( create_persistence_heat_maps_weighted_by_arctan_of_their_persistence create_persistence_heat_maps_weighted_by_arctan_of_their_persistence.cpp )
target_link_libraries(create_persistence_heat_maps_weighted_by_arctan_of_their_persistence ${Boost_SYSTEM_LIBRARY})
add_test ( create_persistence_heat_maps_weighted_by_arctan_of_their_persistence ${CMAKE_CURRENT_BINARY_DIR}/create_persistence_heat_maps_weighted_by_arctan_of_their_persistence 10 -1 -1 4 -1 ../../example/simple_diagram.txt )

add_executable ( average_persistence_heat_maps average_persistence_heat_maps.cpp )
target_link_libraries(average_persistence_heat_maps ${Boost_SYSTEM_LIBRARY})
add_test ( average_persistence_heat_maps ${CMAKE_CURRENT_BINARY_DIR}/average_persistence_heat_maps simple_diagram.txt.mps simple_diagram.txt.mps )

add_executable ( plot_persistence_heat_map plot_persistence_heat_map.cpp )
target_link_libraries(plot_persistence_heat_map ${Boost_SYSTEM_LIBRARY})
add_test ( plot_persistence_heat_map ${CMAKE_CURRENT_BINARY_DIR}/plot_persistence_heat_map simple_diagram.txt.mps )

add_executable ( compute_distance_of_persistence_heat_maps compute_distance_of_persistence_heat_maps.cpp )
target_link_libraries(compute_distance_of_persistence_heat_maps ${Boost_SYSTEM_LIBRARY})
add_test ( compute_distance_of_persistence_heat_maps ${CMAKE_CURRENT_BINARY_DIR}/compute_distance_of_persistence_heat_maps 1 simple_diagram.txt.mps simple_diagram.txt.mps )

add_executable ( compute_scalar_product_of_persistence_heat_maps compute_scalar_product_of_persistence_heat_maps.cpp )
target_link_libraries(compute_scalar_product_of_persistence_heat_maps ${Boost_SYSTEM_LIBRARY})
add_test ( compute_scalar_product_of_persistence_heat_maps ${CMAKE_CURRENT_BINARY_DIR}/compute_scalar_product_of_persistence_heat_maps simple_diagram.txt.mps simple_diagram.txt.mps )

