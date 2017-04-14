cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)

#persitence diagrams
add_executable ( plot_persistence_intervals plot_persistence_intervals.cpp )
target_link_libraries( plot_persistence_intervals ${Boost_SYSTEM_LIBRARY})
add_test ( plot_persistence_intervals ${CMAKE_CURRENT_BINARY_DIR}/plot_persistence_intervals ../../example/simple_diagram.txt )

add_executable ( compute_birth_death_range_in_persistence_diagram compute_birth_death_range_in_persistence_diagram.cpp )
target_link_libraries( compute_birth_death_range_in_persistence_diagram ${Boost_SYSTEM_LIBRARY})
add_test ( compute_birth_death_range_in_persistence_diagram ${CMAKE_CURRENT_BINARY_DIR}/compute_birth_death_range_in_persistence_diagram -1 ../../example/simple_diagram.txt )

add_executable ( compute_number_of_dominant_intervals compute_number_of_dominant_intervals.cpp )
target_link_libraries( compute_number_of_dominant_intervals ${Boost_SYSTEM_LIBRARY})
add_test ( compute_number_of_dominant_intervals ${CMAKE_CURRENT_BINARY_DIR}/compute_number_of_dominant_intervals ../../example/simple_diagram.txt -1 2 )

add_executable ( plot_histogram_of_intervals_lengths plot_histogram_of_intervals_lengths.cpp )
target_link_libraries( plot_histogram_of_intervals_lengths ${Boost_SYSTEM_LIBRARY})
add_test ( plot_histogram_of_intervals_lengths ${CMAKE_CURRENT_BINARY_DIR}/plot_histogram_of_intervals_lengths ../../example/simple_diagram.txt 2 )

add_executable ( plot_persistence_Betti_numbers plot_persistence_Betti_numbers.cpp )
target_link_libraries( plot_persistence_Betti_numbers ${Boost_SYSTEM_LIBRARY})
add_test ( plot_persistence_Betti_numbers ${CMAKE_CURRENT_BINARY_DIR}/plot_persistence_Betti_numbers ../../example/simple_diagram.txt )

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
add_executable ( compute_bottleneck_distance compute_bottleneck_distance.cpp )
target_link_libraries( compute_bottleneck_distance ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
target_link_libraries(compute_bottleneck_distance ${TBB_LIBRARIES})
endif(TBB_FOUND)
add_test ( compute_bottleneck_distance ${CMAKE_CURRENT_BINARY_DIR}/compute_bottleneck_distance ../../example/simple_diagram.txt ../../example/simple_diagram.txt )
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
