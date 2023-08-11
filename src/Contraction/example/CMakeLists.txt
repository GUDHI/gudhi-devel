project(Contraction_examples)

if(TARGET CGAL::CGAL)
  add_executable(RipsContraction Rips_contraction.cpp)
  add_executable(GarlandHeckbert Garland_heckbert.cpp)

  target_link_libraries(RipsContraction CGAL::CGAL)
  target_link_libraries(GarlandHeckbert CGAL::CGAL)

  add_test(NAME Contraction_example_tore3D_0.2 COMMAND $<TARGET_FILE:RipsContraction>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "0.2")

endif(TARGET CGAL::CGAL)
