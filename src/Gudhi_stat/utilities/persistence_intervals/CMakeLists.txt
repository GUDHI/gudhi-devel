cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)

#persitence diagrams
add_executable ( plot_persistence_intervals plot_persistence_intervals.cpp )
target_link_libraries( plot_persistence_intervals ${Boost_SYSTEM_LIBRARY})

add_executable ( compute_birth_death_range_in_persistence_diagram compute_birth_death_range_in_persistence_diagram.cpp )
target_link_libraries( compute_birth_death_range_in_persistence_diagram ${Boost_SYSTEM_LIBRARY})

add_executable ( compute_number_of_dominant_intervals compute_number_of_dominant_intervals.cpp )
target_link_libraries( compute_number_of_dominant_intervals ${Boost_SYSTEM_LIBRARY})

add_executable ( plot_histogram_of_intervals_lengths plot_histogram_of_intervals_lengths.cpp )
target_link_libraries( plot_histogram_of_intervals_lengths ${Boost_SYSTEM_LIBRARY})

add_executable ( plot_persistence_Betti_numbers plot_persistence_Betti_numbers.cpp )
target_link_libraries( plot_persistence_Betti_numbers ${Boost_SYSTEM_LIBRARY})

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
add_executable ( compute_bottleneck_distance compute_bottleneck_distance.cpp )
target_link_libraries( compute_bottleneck_distance ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
target_link_libraries(compute_bottleneck_distance ${TBB_LIBRARIES})
endif(TBB_FOUND)

endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
