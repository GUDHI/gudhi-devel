project(Common_examples)

add_executable ( vector_double_off_reader example_vector_double_points_off_reader.cpp )
target_link_libraries(vector_double_off_reader ${CGAL_LIBRARY})
file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
add_test(NAME Common_example_vector_double_off_reader COMMAND $<TARGET_FILE:vector_double_off_reader>
    "alphacomplexdoc.off")

if (DIFF_PATH)
  # Do not forget to copy test results files in current binary dir
  file(COPY "vectordoubleoffreader_result.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  add_test(Common_example_vector_double_off_reader_diff_files ${DIFF_PATH}
      ${CMAKE_CURRENT_BINARY_DIR}/vectordoubleoffreader_result.txt ${CMAKE_CURRENT_BINARY_DIR}/alphacomplexdoc.off.txt)
  set_tests_properties(Common_example_vector_double_off_reader_diff_files PROPERTIES DEPENDS Common_example_vector_double_off_reader) 
endif()

if(NOT CGAL_VERSION VERSION_LESS 4.11.0)
  add_executable ( cgal_3D_off_reader example_CGAL_3D_points_off_reader.cpp )
  target_link_libraries(cgal_3D_off_reader ${CGAL_LIBRARY})
  add_test(NAME Common_example_vector_cgal_3D_off_reader COMMAND $<TARGET_FILE:cgal_3D_off_reader>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off")
endif()

# requires CGAL and Eigen3
if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  add_executable ( cgal_off_reader example_CGAL_points_off_reader.cpp )
  target_link_libraries(cgal_off_reader ${CGAL_LIBRARY})
  add_test(NAME Common_example_vector_cgal_off_reader COMMAND $<TARGET_FILE:cgal_off_reader>
      "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off")
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
