project(Contraction_examples)

if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
  add_executable(RipsContraction Rips_contraction.cpp)

  add_executable(GarlandHeckbert Garland_heckbert.cpp)
  target_link_libraries(GarlandHeckbert ${Boost_TIMER_LIBRARY})

  add_test(NAME Contraction_example_tore3D_0.2 COMMAND $<TARGET_FILE:RipsContraction>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "0.2")
  # TODO(DS) : These tests are too long under Windows
  #add_test(NAME Contraction_example_sphere_0.2 COMMAND $<TARGET_FILE:RipsContraction>
  #    "${CMAKE_SOURCE_DIR}/data/points/sphere3D_2646.off" "0.2")
  #add_test(NAME Contraction_example_SO3_0.3 COMMAND $<TARGET_FILE:RipsContraction>
  #    "${CMAKE_SOURCE_DIR}/data/points/SO3_10000.off" "0.3")

endif (NOT CGAL_VERSION VERSION_LESS 4.11.0)
