project(Persistence_representations_intervals_utilities)


add_executable ( plot_histogram_of_intervals_lengths plot_histogram_of_intervals_lengths.cpp )

add_test(NAME plot_histogram_of_intervals_lengths COMMAND $<TARGET_FILE:plot_histogram_of_intervals_lengths>
    "${CMAKE_CURRENT_BINARY_DIR}/../first.pers" "-1")

install(TARGETS plot_histogram_of_intervals_lengths DESTINATION bin)

add_persistence_representation_plot_utility(plot_persistence_intervals "")
add_persistence_representation_plot_utility(plot_persistence_Betti_numbers "")

add_persistence_representation_creation_utility(compute_birth_death_range_in_persistence_diagram "-1")


add_executable ( compute_number_of_dominant_intervals compute_number_of_dominant_intervals.cpp )
add_test(NAME Persistence_representation_utilities_compute_number_of_dominant_intervals
    COMMAND $<TARGET_FILE:compute_number_of_dominant_intervals>
    "${CMAKE_CURRENT_BINARY_DIR}/../first.pers" "-1" "2")

install(TARGETS compute_number_of_dominant_intervals DESTINATION bin)


if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  add_executable ( compute_bottleneck_distance compute_bottleneck_distance.cpp )
  if (TBB_FOUND)
  target_link_libraries(compute_bottleneck_distance ${TBB_LIBRARIES})
  endif(TBB_FOUND)
  add_test(NAME Persistence_representation_utilities_compute_bottleneck_distance
      COMMAND $<TARGET_FILE:compute_bottleneck_distance>
      "-1"
      "${CMAKE_CURRENT_BINARY_DIR}/../first.pers"
      "${CMAKE_CURRENT_BINARY_DIR}/../second.pers")

  install(TARGETS compute_bottleneck_distance DESTINATION bin)
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
