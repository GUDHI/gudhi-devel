project(Alpha_complex_examples)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  add_executable ( Alpha_complex_example_from_points Alpha_complex_from_points.cpp )
  target_link_libraries(Alpha_complex_example_from_points  ${CGAL_LIBRARY})
  add_executable ( Alpha_complex_example_from_off Alpha_complex_from_off.cpp )
  target_link_libraries(Alpha_complex_example_from_off ${CGAL_LIBRARY})
  add_executable ( Alpha_complex_example_fast_from_off Fast_alpha_complex_from_off.cpp )
  target_link_libraries(Alpha_complex_example_fast_from_off ${CGAL_LIBRARY})
  if(TARGET TBB::tbb)
    target_link_libraries(Alpha_complex_example_from_points TBB::tbb)
    target_link_libraries(Alpha_complex_example_from_off TBB::tbb)
    target_link_libraries(Alpha_complex_example_fast_from_off TBB::tbb)
  endif()

  add_test(NAME Alpha_complex_example_from_points COMMAND $<TARGET_FILE:Alpha_complex_example_from_points>)

  add_test(NAME Alpha_complex_example_from_off_60 COMMAND $<TARGET_FILE:Alpha_complex_example_from_off>
      "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" "60.0" "${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_result_60.txt")
  add_test(NAME Alpha_complex_example_from_off_32 COMMAND $<TARGET_FILE:Alpha_complex_example_from_off>
      "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" "32.0" "${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_result_32.txt")

  add_test(NAME Alpha_complex_example_fast_from_off_60 COMMAND $<TARGET_FILE:Alpha_complex_example_fast_from_off>
      "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" "60.0" "${CMAKE_CURRENT_BINARY_DIR}/fastalphaoffreader_result_60.txt")
  add_test(NAME Alpha_complex_example_fast_from_off_32 COMMAND $<TARGET_FILE:Alpha_complex_example_fast_from_off>
      "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" "32.0" "${CMAKE_CURRENT_BINARY_DIR}/fastalphaoffreader_result_32.txt")

  if (DIFF_PATH)
    # Do not forget to copy test results files in current binary dir
    file(COPY "alphaoffreader_for_doc_32.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
    file(COPY "alphaoffreader_for_doc_60.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

    add_test(Alpha_complex_example_from_off_60_diff_files ${DIFF_PATH}
        ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_result_60.txt ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_for_doc_60.txt)
    set_tests_properties(Alpha_complex_example_from_off_60_diff_files PROPERTIES DEPENDS Alpha_complex_example_from_off_60) 
    add_test(Alpha_complex_example_from_off_32_diff_files ${DIFF_PATH}
        ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_result_32.txt ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_for_doc_32.txt)
    set_tests_properties(Alpha_complex_example_from_off_32_diff_files PROPERTIES DEPENDS Alpha_complex_example_from_off_32) 

    add_test(Alpha_complex_example_fast_from_off_60_diff_files ${DIFF_PATH}
        ${CMAKE_CURRENT_BINARY_DIR}/fastalphaoffreader_result_60.txt ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_for_doc_60.txt)
    set_tests_properties(Alpha_complex_example_fast_from_off_60_diff_files PROPERTIES DEPENDS Alpha_complex_example_fast_from_off_60) 
    add_test(Alpha_complex_example_fast_from_off_32_diff_files ${DIFF_PATH}
        ${CMAKE_CURRENT_BINARY_DIR}/fastalphaoffreader_result_32.txt ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_for_doc_32.txt)
    set_tests_properties(Alpha_complex_example_fast_from_off_32_diff_files PROPERTIES DEPENDS Alpha_complex_example_fast_from_off_32) 
  endif()
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)

if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
  add_executable ( Alpha_complex_example_weighted_3d_from_points Weighted_alpha_complex_3d_from_points.cpp )
  target_link_libraries(Alpha_complex_example_weighted_3d_from_points  ${CGAL_LIBRARY})
  if(TARGET TBB::tbb)
    target_link_libraries(Alpha_complex_example_weighted_3d_from_points TBB::tbb)
  endif()
  add_test(NAME Alpha_complex_example_weighted_3d_from_points
          COMMAND $<TARGET_FILE:Alpha_complex_example_weighted_3d_from_points>)

  add_executable ( Alpha_complex_example_3d_from_points Alpha_complex_3d_from_points.cpp )
  target_link_libraries(Alpha_complex_example_3d_from_points  ${CGAL_LIBRARY})
  if(TARGET TBB::tbb)
    target_link_libraries(Alpha_complex_example_3d_from_points TBB::tbb)
  endif()
  add_test(NAME Alpha_complex_example_3d_from_points
          COMMAND $<TARGET_FILE:Alpha_complex_example_3d_from_points>)

endif(NOT CGAL_VERSION VERSION_LESS 4.11.0)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 5.1.0)
  add_executable ( Weighted_alpha_complex_example_from_points Weighted_alpha_complex_from_points.cpp )
  target_link_libraries(Weighted_alpha_complex_example_from_points  ${CGAL_LIBRARY})
  if(TARGET TBB::tbb)
    target_link_libraries(Weighted_alpha_complex_example_from_points TBB::tbb)
  endif()
  add_test(NAME Weighted_alpha_complex_example_from_points COMMAND $<TARGET_FILE:Weighted_alpha_complex_example_from_points>)
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 5.1.0)
