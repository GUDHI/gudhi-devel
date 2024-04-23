add_executable ( vector_double_off_reader example_vector_double_points_off_reader.cpp )
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

add_executable_with_targets(cgal_3D_off_reader example_CGAL_3D_points_off_reader.cpp CGAL::CGAL TBB::tbb)

if(TARGET CGAL::CGAL)
  add_test(NAME Common_example_vector_cgal_3D_off_reader COMMAND $<TARGET_FILE:cgal_3D_off_reader>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off")
endif()

add_executable_with_targets(cgal_off_reader example_CGAL_points_off_reader.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  add_test(NAME Common_example_vector_cgal_off_reader COMMAND $<TARGET_FILE:cgal_off_reader>
      "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off")
endif()
