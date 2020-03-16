project(Persistence_representations_heat_maps_utilities)

add_persistence_representation_creation_utility(create_pssk  "10" "-1" "-1" "4" "-1")
add_persistence_representation_creation_utility(create_p_h_m_weighted_by_arctan_of_their_persistence "10" "-1" "-1" "4" "-1")

add_persistence_representation_creation_utility(create_p_h_m_weighted_by_distance_from_diagonal "10" "-1" "-1" "4" "-1")
# Tests output the same file
set_tests_properties(Persistence_representation_utilities_create_p_h_m_weighted_by_distance_from_diagonal PROPERTIES DEPENDS
                     Persistence_representation_utilities_create_p_h_m_weighted_by_arctan_of_their_persistence)

add_persistence_representation_creation_utility(create_p_h_m_weighted_by_squared_diag_distance "10" "-1" "-1" "4" "-1")
# Tests output the same file
set_tests_properties(Persistence_representation_utilities_create_p_h_m_weighted_by_squared_diag_distance PROPERTIES DEPENDS
                     Persistence_representation_utilities_create_p_h_m_weighted_by_distance_from_diagonal)

# Need to set grid min and max for further average, distance and scalar_product
add_persistence_representation_creation_utility(create_persistence_heat_maps "10" "0" "35" "10" "-1")
set_tests_properties(Persistence_representation_utilities_create_persistence_heat_maps PROPERTIES DEPENDS
                     Persistence_representation_utilities_create_p_h_m_weighted_by_squared_diag_distance)

add_persistence_representation_plot_utility(create_persistence_heat_maps plot_persistence_heat_map ".mps")

add_persistence_representation_function_utility(create_persistence_heat_maps average_persistence_heat_maps ".mps")
add_persistence_representation_function_utility(create_persistence_heat_maps compute_distance_of_persistence_heat_maps ".mps" "1")
add_persistence_representation_function_utility(create_persistence_heat_maps compute_scalar_product_of_persistence_heat_maps ".mps")
