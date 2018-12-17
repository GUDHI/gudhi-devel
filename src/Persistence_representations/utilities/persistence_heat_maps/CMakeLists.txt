project(Persistence_representations_heat_maps_utilities)

add_persistence_representation_creation_utility(create_pssk  "10" "-1" "-1" "4" "-1")
add_persistence_representation_creation_utility(create_p_h_m_weighted_by_arctan_of_their_persistence "10" "-1" "-1" "4" "-1")
add_persistence_representation_creation_utility(create_p_h_m_weighted_by_distance_from_diagonal "10" "-1" "-1" "4" "-1")
add_persistence_representation_creation_utility(create_p_h_m_weighted_by_squared_diag_distance "10" "-1" "-1" "4" "-1")
# Need to set grid min and max for further average, distance and scalar_product
add_persistence_representation_creation_utility(create_persistence_heat_maps "10" "0" "35" "10" "-1")

add_persistence_representation_plot_utility(plot_persistence_heat_map ".mps")

add_persistence_representation_function_utility(average_persistence_heat_maps ".mps")
add_persistence_representation_function_utility(compute_distance_of_persistence_heat_maps ".mps" "1")
add_persistence_representation_function_utility(compute_scalar_product_of_persistence_heat_maps ".mps")
