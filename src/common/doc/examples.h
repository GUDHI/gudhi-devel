// List of GUDHI examples and utils - Doxygen needs at least a file tag to analyse comments
// Generated from scripts/cpp_examples_for_doxygen.py
/*! @file Examples
 * \section Skeleton_blocker_example_section Skeleton_blocker
 * @example Skeleton_blocker_iteration.cpp
 * @example Skeleton_blocker_from_simplices.cpp
 * @example Skeleton_blocker_link.cpp
 * \section Alpha_complex_example_section Alpha_complex
 * @example alpha_complex_persistence.cpp
 * @example alpha_complex_3d_persistence.cpp
 * @example Alpha_complex_from_points.cpp
 * @example Alpha_complex_3d_from_points.cpp
 * @example Weighted_alpha_complex_3d_from_points.cpp
 * @example Fast_alpha_complex_from_off.cpp
 * @example Alpha_complex_from_off.cpp
 * @example Weighted_alpha_complex_from_points.cpp
 * \section Simplex_tree_example_section Simplex_tree
 * @example cech_complex_cgal_mini_sphere_3d.cpp
 * @example mini_simplex_tree.cpp
 * @example example_alpha_shapes_3_simplex_tree_from_off_file.cpp
 * @example graph_expansion_with_blocker.cpp
 * @example simple_simplex_tree.cpp
 * @example simplex_tree_from_cliques_of_graph.cpp
 * \section Collapse_example_section Collapse
 * @example point_cloud_edge_collapse_rips_persistence.cpp
 * @example distance_matrix_edge_collapse_rips_persistence.cpp
 * @example edge_collapse_basic_example.cpp
 * @example edge_collapse_conserve_persistence.cpp
 * \section Persistent_cohomology_example_section Persistent_cohomology
 * @example rips_persistence_via_boundary_matrix.cpp
 * @example rips_persistence_step_by_step.cpp
 * @example plain_homology.cpp
 * @example custom_persistence_sort.cpp
 * @example persistence_from_simple_simplex_tree.cpp
 * @example persistence_from_file.cpp
 * @example rips_multifield_persistence.cpp
 * \section Nerve_GIC_example_section Nerve_GIC
 * @example VoronoiGIC.cpp
 * @example Nerve.cpp
 * @example FuncGIC.cpp
 * @example CoordGIC.cpp
 * \section Rips_complex_example_section Rips_complex
 * @example rips_persistence.cpp
 * @example sparse_rips_persistence.cpp
 * @example rips_correlation_matrix_persistence.cpp
 * @example rips_distance_matrix_persistence.cpp
 * @example example_one_skeleton_rips_from_correlation_matrix.cpp
 * @example example_rips_complex_from_csv_distance_matrix_file.cpp
 * @example example_sparse_rips.cpp
 * @example example_one_skeleton_rips_from_distance_matrix.cpp
 * @example example_rips_complex_from_off_file.cpp
 * @example example_one_skeleton_rips_from_points.cpp
 * \section Persistence_representations_example_section Persistence_representations
 * @example persistence_landscapes_on_grid/average_landscapes_on_grid.cpp
 * @example persistence_landscapes_on_grid/create_landscapes_on_grid.cpp
 * @example persistence_landscapes_on_grid/compute_distance_of_landscapes_on_grid.cpp
 * @example persistence_landscapes_on_grid/plot_landscapes_on_grid.cpp
 * @example persistence_landscapes_on_grid/compute_scalar_product_of_landscapes_on_grid.cpp
 * @example persistence_intervals/compute_birth_death_range_in_persistence_diagram.cpp
 * @example persistence_intervals/plot_persistence_intervals.cpp
 * @example persistence_intervals/plot_histogram_of_intervals_lengths.cpp
 * @example persistence_intervals/compute_bottleneck_distance.cpp
 * @example persistence_intervals/plot_persistence_Betti_numbers.cpp
 * @example persistence_intervals/compute_number_of_dominant_intervals.cpp
 * @example persistence_heat_maps/average_persistence_heat_maps.cpp
 * @example persistence_heat_maps/create_p_h_m_weighted_by_squared_diag_distance.cpp
 * @example persistence_heat_maps/create_pssk.cpp
 * @example persistence_heat_maps/create_p_h_m_weighted_by_distance_from_diagonal.cpp
 * @example persistence_heat_maps/compute_distance_of_persistence_heat_maps.cpp
 * @example persistence_heat_maps/create_p_h_m_weighted_by_arctan_of_their_persistence.cpp
 * @example persistence_heat_maps/plot_persistence_heat_map.cpp
 * @example persistence_heat_maps/compute_scalar_product_of_persistence_heat_maps.cpp
 * @example persistence_heat_maps/create_persistence_heat_maps.cpp
 * @example persistence_vectors/average_persistence_vectors.cpp
 * @example persistence_vectors/plot_persistence_vectors.cpp
 * @example persistence_vectors/create_persistence_vectors.cpp
 * @example persistence_vectors/compute_scalar_product_of_persistence_vectors.cpp
 * @example persistence_vectors/compute_distance_of_persistence_vectors.cpp
 * @example persistence_landscapes/create_landscapes.cpp
 * @example persistence_landscapes/average_landscapes.cpp
 * @example persistence_landscapes/compute_distance_of_landscapes.cpp
 * @example persistence_landscapes/plot_landscapes.cpp
 * @example persistence_landscapes/compute_scalar_product_of_landscapes.cpp
 * @example persistence_heat_maps.cpp
 * @example sliced_wasserstein.cpp
 * @example persistence_intervals.cpp
 * @example persistence_vectors.cpp
 * @example persistence_landscape.cpp
 * @example persistence_landscape_on_grid.cpp
 * \section common_example_section common
 * @example off_file_from_shape_generator.cpp
 * @example example_CGAL_points_off_reader.cpp
 * @example example_vector_double_points_off_reader.cpp
 * @example example_CGAL_3D_points_off_reader.cpp
 * \section Subsampling_example_section Subsampling
 * @example example_pick_n_random_points.cpp
 * @example example_choose_n_farthest_points.cpp
 * @example example_sparsify_point_set.cpp
 * @example example_custom_distance.cpp
 * \section Contraction_example_section Contraction
 * @example Rips_contraction.cpp
 * @example Garland_heckbert.cpp
 * \section Bottleneck_distance_example_section Bottleneck_distance
 * @example bottleneck_distance.cpp
 * @example bottleneck_basic_example.cpp
 * @example alpha_rips_persistence_bottleneck_distance.cpp
 * \section Tangential_complex_example_section Tangential_complex
 * @example example_basic.cpp
 * @example example_with_perturb.cpp
 * \section Cech_complex_example_section Cech_complex
 * @example cech_persistence.cpp
 * @example cech_complex_example_from_points.cpp
 * @example cech_complex_step_by_step.cpp
 * \section Spatial_searching_example_section Spatial_searching
 * @example example_spatial_searching.cpp
 * \section Bitmap_cubical_complex_example_section Bitmap_cubical_complex
 * @example periodic_cubical_complex_persistence.cpp
 * @example cubical_complex_persistence.cpp
 * @example Random_bitmap_cubical_complex.cpp
 * \section Toplex_map_example_section Toplex_map
 * @example simple_toplex_map.cpp
 * \section Witness_complex_example_section Witness_complex
 * @example weak_witness_persistence.cpp
 * @example strong_witness_persistence.cpp
 * @example example_witness_complex_sphere.cpp
 * @example example_strong_witness_complex_off.cpp
 * @example example_nearest_landmark_table.cpp
 * @example example_witness_complex_off.cpp
 */
 
