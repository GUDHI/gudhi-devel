cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)



#persitence landscapes
add_executable ( create_landscapes create_landscapes.cpp )
target_link_libraries(create_landscapes ${Boost_SYSTEM_LIBRARY})
add_test ( create_landscapes ${CMAKE_CURRENT_BINARY_DIR}/create_landscapes -1 ../../example/simple_diagram.txt )

add_executable ( average_landscapes average_landscapes.cpp )
target_link_libraries(average_landscapes ${Boost_SYSTEM_LIBRARY})
add_test ( average_landscapes ${CMAKE_CURRENT_BINARY_DIR}/average_landscapes simple_diagram.txt.land simple_diagram.txt.land simple_diagram.txt.land )

add_executable ( plot_landscapes plot_landscapes.cpp )
target_link_libraries(plot_landscapes ${Boost_SYSTEM_LIBRARY})
add_test ( plot_landscapes ${CMAKE_CURRENT_BINARY_DIR}/plot_landscapes simple_diagram.txt.land )

add_executable ( compute_distance_of_landscapes compute_distance_of_landscapes.cpp )
target_link_libraries(compute_distance_of_landscapes ${Boost_SYSTEM_LIBRARY})
add_test ( compute_distance_of_landscapes ${CMAKE_CURRENT_BINARY_DIR}/compute_distance_of_landscapes 1 simple_diagram.txt.land simple_diagram.txt.land )

add_executable ( compute_scalar_product_of_landscapes compute_scalar_product_of_landscapes.cpp )
target_link_libraries(compute_scalar_product_of_landscapes ${Boost_SYSTEM_LIBRARY})
add_test ( compute_scalar_product_of_landscapes ${CMAKE_CURRENT_BINARY_DIR}/compute_scalar_product_of_landscapes simple_diagram.txt.land simple_diagram.txt.land )

