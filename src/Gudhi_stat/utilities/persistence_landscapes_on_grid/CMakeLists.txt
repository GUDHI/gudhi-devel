cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)




#persitence landscapes on a grid
add_executable ( create_landscapes_on_grid create_landscapes_on_grid.cpp )
target_link_libraries(create_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})
add_test ( create_landscapes_on_grid ${CMAKE_CURRENT_BINARY_DIR}/create_landscapes_on_grid 100 -1 -1 -1 "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt" )

add_executable ( average_landscapes_on_grid average_landscapes_on_grid.cpp )
target_link_libraries(average_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})
add_test ( average_landscapes_on_grid ${CMAKE_CURRENT_BINARY_DIR}/average_landscapes_on_grid "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.g_land" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.g_land" )

add_executable ( plot_landscapes_on_grid plot_landscapes_on_grid.cpp )
target_link_libraries(plot_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})
add_test ( plot_landscapes_on_grid ${CMAKE_CURRENT_BINARY_DIR}/plot_landscapes_on_grid "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.g_land" )

add_executable ( compute_distance_of_landscapes_on_grid compute_distance_of_landscapes_on_grid.cpp )
target_link_libraries(compute_distance_of_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})
add_test ( compute_distance_of_landscapes_on_grids ${CMAKE_CURRENT_BINARY_DIR}/compute_distance_of_landscapes_on_grid "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.g_land" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.g_land" )

add_executable ( compute_scalar_product_of_landscapes_on_grid compute_scalar_product_of_landscapes_on_grid.cpp )
target_link_libraries(compute_scalar_product_of_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})
add_test ( compute_scalar_product_of_landscapes_on_grid ${CMAKE_CURRENT_BINARY_DIR}/compute_scalar_product_of_landscapes_on_grid "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.g_land" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt.g_land" )	

