project(Witness_complex_utilities)

# CGAL and Eigen3 are required for Euclidean version of Witness
if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  if(TARGET Boost::program_options)
    add_executable ( Witness_complex_strong_witness_persistence strong_witness_persistence.cpp )
    target_link_libraries(Witness_complex_strong_witness_persistence Boost::program_options)

    add_executable ( Witness_complex_weak_witness_persistence weak_witness_persistence.cpp )
    target_link_libraries(Witness_complex_weak_witness_persistence Boost::program_options)

    if(TARGET TBB::tbb)
      target_link_libraries(Witness_complex_strong_witness_persistence TBB::tbb)
      target_link_libraries(Witness_complex_weak_witness_persistence TBB::tbb)
    endif()

    add_test(NAME Witness_complex_strong_test_torus_persistence
        COMMAND $<TARGET_FILE:Witness_complex_strong_witness_persistence>
        "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-l" "20" "-a" "0.5")
    add_test(NAME Witness_complex_weak_test_torus_persistence
        COMMAND $<TARGET_FILE:Witness_complex_weak_witness_persistence>
        "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-l" "20" "-a" "0.5")

    install(TARGETS Witness_complex_strong_witness_persistence DESTINATION bin)
    install(TARGETS Witness_complex_weak_witness_persistence DESTINATION bin)
  endif()
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
