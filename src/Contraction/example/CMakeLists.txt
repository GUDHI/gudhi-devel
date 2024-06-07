add_executable_with_targets(RipsContraction Rips_contraction.cpp CGAL::CGAL TBB::tbb)
add_executable_with_targets(GarlandHeckbert Garland_heckbert.cpp CGAL::CGAL TBB::tbb)

if(TARGET CGAL::CGAL)
  add_test(NAME Contraction_example_tore3D_0.2 COMMAND $<TARGET_FILE:RipsContraction>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "0.2")
endif()
