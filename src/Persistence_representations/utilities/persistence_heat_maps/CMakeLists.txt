cmake_minimum_required(VERSION 2.6)
project(Persistence_representations_heat_maps_utilities)

file(COPY "${CMAKE_SOURCE_DIR}/data/persistence_diagram/simple_diagram.txt" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

add_executable ( create_persistence_heat_maps create_persistence_heat_maps.cpp )
target_link_libraries(create_persistence_heat_maps ${Boost_SYSTEM_LIBRARY})

add_test(NAME create_persistence_heat_maps COMMAND $<TARGET_FILE:create_persistence_heat_maps>
    "10" "-1" "-1" "4" "-1" "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")


add_executable ( create_pssk create_pssk.cpp )
target_link_libraries(create_pssk ${Boost_SYSTEM_LIBRARY})

add_test(NAME create_pssk COMMAND $<TARGET_FILE:create_pssk>
    "10" "-1" "-1" "4" "-1" "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( create_p_h_m_weighted_by_distance_from_diagonal create_p_h_m_weighted_by_distance_from_diagonal.cpp  )
target_link_libraries(create_p_h_m_weighted_by_distance_from_diagonal ${Boost_SYSTEM_LIBRARY})

add_test(NAME create_p_h_m_weighted_by_distance_from_diagonal COMMAND $<TARGET_FILE:create_p_h_m_weighted_by_distance_from_diagonal>
    "10" "-1" "-1" "4" "-1" "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( create_p_h_m_weighted_by_squared_diag_distance create_p_h_m_weighted_by_squared_diag_distance.cpp )
target_link_libraries(create_p_h_m_weighted_by_squared_diag_distance ${Boost_SYSTEM_LIBRARY})

add_test(NAME create_p_h_m_weighted_by_squared_diag_distance COMMAND $<TARGET_FILE:create_p_h_m_weighted_by_squared_diag_distance>
    "10" "-1" "-1" "4" "-1" "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( create_p_h_m_weighted_by_arctan_of_their_persistence create_p_h_m_weighted_by_arctan_of_their_persistence.cpp )
target_link_libraries(create_p_h_m_weighted_by_arctan_of_their_persistence ${Boost_SYSTEM_LIBRARY})

add_test(NAME create_p_h_m_weighted_by_arctan_of_their_persistence COMMAND $<TARGET_FILE:create_p_h_m_weighted_by_arctan_of_their_persistence>
    "10" "-1" "-1" "4" "-1" "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( average_persistence_heat_maps average_persistence_heat_maps.cpp )
target_link_libraries(average_persistence_heat_maps ${Boost_SYSTEM_LIBRARY})

add_test(NAME average_persistence_heat_maps COMMAND $<TARGET_FILE:average_persistence_heat_maps>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( plot_persistence_heat_map plot_persistence_heat_map.cpp )
target_link_libraries(plot_persistence_heat_map ${Boost_SYSTEM_LIBRARY})

add_test(NAME plot_persistence_heat_map COMMAND $<TARGET_FILE:plot_persistence_heat_map>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( compute_distance_of_persistence_heat_maps compute_distance_of_persistence_heat_maps.cpp )
target_link_libraries(compute_distance_of_persistence_heat_maps ${Boost_SYSTEM_LIBRARY})

add_test(NAME compute_distance_of_persistence_heat_maps COMMAND $<TARGET_FILE:compute_distance_of_persistence_heat_maps>
    "1"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")

add_executable ( compute_scalar_product_of_persistence_heat_maps compute_scalar_product_of_persistence_heat_maps.cpp )
target_link_libraries(compute_scalar_product_of_persistence_heat_maps ${Boost_SYSTEM_LIBRARY})

add_test(NAME compute_scalar_product_of_persistence_heat_maps COMMAND $<TARGET_FILE:compute_scalar_product_of_persistence_heat_maps>
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt"
    "${CMAKE_CURRENT_BINARY_DIR}/simple_diagram.txt")
