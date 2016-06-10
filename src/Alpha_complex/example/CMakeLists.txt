cmake_minimum_required(VERSION 2.6)
project(Alpha_shapes_examples)

# need CGAL 4.7
# cmake -DCGAL_DIR=~/workspace/CGAL-4.7-Ic-41 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )

      add_executable ( alphapoints Alpha_complex_from_points.cpp )
      target_link_libraries(alphapoints  ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${CGAL_LIBRARY})
      add_executable ( alphaoffreader Alpha_complex_from_off.cpp )
      target_link_libraries(alphaoffreader  ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${CGAL_LIBRARY})
      if (TBB_FOUND)
        target_link_libraries(alphapoints ${TBB_LIBRARIES})
        target_link_libraries(alphaoffreader ${TBB_LIBRARIES})
      endif()

      add_test(alphapoints ${CMAKE_CURRENT_BINARY_DIR}/alphapoints)
      # Do not forget to copy test files in current binary dir
      file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
      add_test(alphaoffreader_doc_60 ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader alphacomplexdoc.off 60.0 ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_result_60.txt)
      add_test(alphaoffreader_doc_32 ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader alphacomplexdoc.off 32.0 ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_result_32.txt)
      if (DIFF_PATH)
        # Do not forget to copy test results files in current binary dir
        file(COPY "alphaoffreader_for_doc_32.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
        file(COPY "alphaoffreader_for_doc_60.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

        add_test(alphaoffreader_doc_60_diff_files ${DIFF_PATH} ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_result_60.txt ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_for_doc_60.txt)
        add_test(alphaoffreader_doc_32_diff_files ${DIFF_PATH} ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_result_32.txt ${CMAKE_CURRENT_BINARY_DIR}/alphaoffreader_for_doc_32.txt)
      endif()
    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()
