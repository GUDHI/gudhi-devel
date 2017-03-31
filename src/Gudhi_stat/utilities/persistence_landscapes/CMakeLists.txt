cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)



#persitence landscapes
add_executable ( average_landscapes average_landscapes.cpp )
target_link_libraries(average_landscapes ${Boost_SYSTEM_LIBRARY})

add_executable ( create_landscapes create_landscapes.cpp )
target_link_libraries(create_landscapes ${Boost_SYSTEM_LIBRARY})

add_executable ( plot_landscapes plot_landscapes.cpp )
target_link_libraries(plot_landscapes ${Boost_SYSTEM_LIBRARY})

add_executable ( compute_distance_of_landscapes compute_distance_of_landscapes.cpp )
target_link_libraries(compute_distance_of_landscapes ${Boost_SYSTEM_LIBRARY})

add_executable ( compute_scalar_product_of_landscapes compute_scalar_product_of_landscapes.cpp )
target_link_libraries(compute_scalar_product_of_landscapes ${Boost_SYSTEM_LIBRARY})

