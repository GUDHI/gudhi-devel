project(Nerve_GIC_examples)

if (NOT CGAL_VERSION VERSION_LESS 4.11.0)

  add_executable ( CoordGIC CoordGIC.cpp )
  add_executable ( FuncGIC FuncGIC.cpp )

  if (TBB_FOUND)
    target_link_libraries(CoordGIC ${TBB_LIBRARIES})
    target_link_libraries(FuncGIC ${TBB_LIBRARIES})
  endif()

  # Copy files for not to pollute sources when testing
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/COIL_database/lucky_cat.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/COIL_database/lucky_cat_PCA1" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  add_test(NAME Nerve_GIC_example_CoordGIC COMMAND $<TARGET_FILE:CoordGIC>
      "${CMAKE_CURRENT_BINARY_DIR}/tore3D_1307.off" "0")

  add_test(NAME Nerve_GIC_example_FuncGIC COMMAND $<TARGET_FILE:FuncGIC>
      "${CMAKE_CURRENT_BINARY_DIR}/lucky_cat.off"
      "${CMAKE_CURRENT_BINARY_DIR}/lucky_cat_PCA1")

  install(TARGETS CoordGIC DESTINATION bin)
  install(TARGETS FuncGIC DESTINATION bin)

endif (NOT CGAL_VERSION VERSION_LESS 4.11.0)
