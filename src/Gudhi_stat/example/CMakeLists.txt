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


##UTILITIES:

#persitence diagrams
add_executable ( utilities/persistence_intervals/plot_persistence_intervals utilities/persistence_intervals/plot_persistence_intervals.cpp )
target_link_libraries( utilities/persistence_intervals/plot_persistence_intervals ${Boost_SYSTEM_LIBRARY})

add_executable ( utilities/persistence_intervals/compute_birth_death_range_in_persistence_diagram utilities/persistence_intervals/compute_birth_death_range_in_persistence_diagram.cpp )
target_link_libraries( utilities/persistence_intervals/compute_birth_death_range_in_persistence_diagram ${Boost_SYSTEM_LIBRARY})

add_executable ( utilities/persistence_intervals/compute_number_of_dominant_intervals utilities/persistence_intervals/compute_number_of_dominant_intervals.cpp )
target_link_libraries( utilities/persistence_intervals/compute_number_of_dominant_intervals ${Boost_SYSTEM_LIBRARY})

add_executable ( utilities/persistence_intervals/plot_histogram_of_intervals_lengths utilities/persistence_intervals/plot_histogram_of_intervals_lengths.cpp )
target_link_libraries( utilities/persistence_intervals/plot_histogram_of_intervals_lengths ${Boost_SYSTEM_LIBRARY})

add_executable ( utilities/persistence_intervals/plot_persistence_Betti_numebrs utilities/persistence_intervals/plot_persistence_Betti_numebrs.cpp )
target_link_libraries( utilities/persistence_intervals/plot_persistence_Betti_numebrs ${Boost_SYSTEM_LIBRARY})

#when we have bottleneck and wasserstein distance computations, add suitable programs here.


#persitence landscapes
add_executable ( utilities/persistence_landscapes/average_landscapes utilities/persistence_landscapes/average_landscapes.cpp )
target_link_libraries(utilities/persistence_landscapes/average_landscapes ${Boost_SYSTEM_LIBRARY})

add_executable ( utilities/persistence_landscapes/create_landscapes utilities/persistence_landscapes/create_landscapes.cpp )
target_link_libraries(utilities/persistence_landscapes/create_landscapes ${Boost_SYSTEM_LIBRARY})

add_executable ( utilities/persistence_landscapes/plot_landscapes utilities/persistence_landscapes/plot_landscapes.cpp )
target_link_libraries(utilities/persistence_landscapes/plot_landscapes ${Boost_SYSTEM_LIBRARY})

#persitence landscapes on a grid

#persitence images

#persitence vectors
