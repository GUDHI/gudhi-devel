cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)




#persitence landscapes on a grid
add_executable ( average_landscapes_on_grid average_landscapes_on_grid.cpp )
target_link_libraries(average_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})

add_executable ( create_landscapes_on_grid create_landscapes_on_grid.cpp )
target_link_libraries(create_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})

add_executable ( plot_landscapes_on_grid plot_landscapes_on_grid.cpp )
target_link_libraries(plot_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})

add_executable ( compute_distance_of_landscapes_on_grid compute_distance_of_landscapes_on_grid.cpp )
target_link_libraries(compute_distance_of_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})

add_executable ( compute_scalar_product_of_landscapes_on_grid compute_scalar_product_of_landscapes_on_grid.cpp )
target_link_libraries(compute_scalar_product_of_landscapes_on_grid ${Boost_SYSTEM_LIBRARY})

