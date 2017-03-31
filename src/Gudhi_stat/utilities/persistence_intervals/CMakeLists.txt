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

add_executable ( plot_persistence_Betti_numebrs plot_persistence_Betti_numebrs.cpp )
target_link_libraries( plot_persistence_Betti_numebrs ${Boost_SYSTEM_LIBRARY})

#when we have bottleneck and wasserstein distance computations, add suitable programs here.
