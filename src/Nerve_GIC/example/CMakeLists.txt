project(Nerve_GIC_examples)

add_executable ( CoordGIC CoordGIC.cpp )
add_executable ( FuncGIC FuncGIC.cpp )

if(TARGET TBB::tbb)
  target_link_libraries(CoordGIC TBB::tbb)
  target_link_libraries(FuncGIC TBB::tbb)
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
