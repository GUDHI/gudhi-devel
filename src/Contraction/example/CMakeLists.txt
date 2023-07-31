project(Contraction_examples)

if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
  add_executable(RipsContraction Rips_contraction.cpp)

  add_executable(GarlandHeckbert Garland_heckbert.cpp)

  add_test(NAME Contraction_example_tore3D_0.2 COMMAND $<TARGET_FILE:RipsContraction>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "0.2")

endif (NOT CGAL_VERSION VERSION_LESS 4.11.0)
