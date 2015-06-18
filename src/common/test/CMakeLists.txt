cmake_minimum_required(VERSION 2.6)
project(GUDHIDelaunayTriangulationOffFileReadWriteUT)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} --coverage")
endif()

# need CGAL 4.6
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )
 
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      
      add_executable ( dtoffrw_UT dtoffrw_unit_test.cpp )
      target_link_libraries(dtoffrw_UT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

      # Unitary tests
      add_test(dtoffrw_UT ${CMAKE_CURRENT_BINARY_DIR}/dtoffrw_UT)

      if (DIFF_PATH)
        add_test(diff_files_UT ${DIFF_PATH} ${CMAKE_CURRENT_BINARY_DIR}/UT.off ${CMAKE_CURRENT_BINARY_DIR}/dtoffrw_alphashapedoc_result.off)
      endif()

    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()

