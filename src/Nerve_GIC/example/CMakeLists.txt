cmake_minimum_required(VERSION 2.6)
project(Nerve_GIC_examples)

if (NOT CGAL_VERSION VERSION_LESS 4.8.1)

  add_executable ( Nerve Nerve.cpp )
  add_executable ( CoordGIC CoordGIC.cpp )
  add_executable ( FuncGIC FuncGIC.cpp )
  add_executable ( VoronoiGIC VoronoiGIC.cpp )

  if (TBB_FOUND)
    target_link_libraries(Nerve ${TBB_LIBRARIES})
    target_link_libraries(CoordGIC ${TBB_LIBRARIES})
    target_link_libraries(FuncGIC ${TBB_LIBRARIES})
    target_link_libraries(VoronoiGIC ${TBB_LIBRARIES})
  endif()

  file(COPY KeplerMapperVisuFromTxtFile.py km.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  add_test(NAME Nerve_GIC_example_nerve COMMAND $<TARGET_FILE:Nerve>
      "${CMAKE_SOURCE_DIR}/data/points/human.off" "2" "10" "0.3")

  add_test(NAME Nerve_GIC_example_VoronoiGIC COMMAND $<TARGET_FILE:VoronoiGIC>
      "${CMAKE_SOURCE_DIR}/data/points/human.off" "100")

  add_test(NAME Nerve_GIC_example_CoordGIC COMMAND $<TARGET_FILE:CoordGIC>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "0")

  add_test(NAME Nerve_GIC_example_FuncGIC COMMAND $<TARGET_FILE:FuncGIC>
      "${CMAKE_SOURCE_DIR}/data/points/COIL_database/lucky_cat.off"
      "${CMAKE_SOURCE_DIR}/data/points/COIL_database/lucky_cat_PCA1")

endif (NOT CGAL_VERSION VERSION_LESS 4.8.1)
