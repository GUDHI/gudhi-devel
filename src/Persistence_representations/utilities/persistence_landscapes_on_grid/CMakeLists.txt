project(Persistence_representations_lanscapes_on_grid_utilities)

# Need to set grid min and max for further average, distance and scalar_product
add_persistence_representation_creation_utility(create_landscapes_on_grid "100" "0" "35" "-1")

add_persistence_representation_plot_utility(plot_landscapes_on_grid ".g_land")

add_persistence_representation_function_utility(average_landscapes_on_grid ".g_land")
add_persistence_representation_function_utility(compute_distance_of_landscapes_on_grid ".g_land" "1")
add_persistence_representation_function_utility(compute_scalar_product_of_landscapes_on_grid ".g_land")
