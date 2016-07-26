cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)

add_executable ( persistence_landscape persistence_landscape.cpp )
target_link_libraries(persistence_landscape ${Boost_SYSTEM_LIBRARY})

add_executable ( average_landscapes average_landscapes.cpp )
target_link_libraries(average_landscapes ${Boost_SYSTEM_LIBRARY})

add_executable ( create_landscapes create_landscapes.cpp )
target_link_libraries(create_landscapes ${Boost_SYSTEM_LIBRARY})

add_executable ( plot_landscapes plot_landscapes.cpp )
target_link_libraries(plot_landscapes ${Boost_SYSTEM_LIBRARY})


add_executable ( persistence_intervals persistence_intervals.cpp )
target_link_libraries(persistence_intervals ${Boost_SYSTEM_LIBRARY})

add_executable ( plot_persistence_intervals plot_persistence_intervals.cpp )
target_link_libraries( plot_persistence_intervals ${Boost_SYSTEM_LIBRARY})

add_executable ( vector_representation vector_representation.cpp )
target_link_libraries(vector_representation ${Boost_SYSTEM_LIBRARY})

