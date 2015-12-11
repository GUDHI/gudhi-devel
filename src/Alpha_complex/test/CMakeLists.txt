cmake_minimum_required(VERSION 2.6)
project(GUDHIAlphaComplexTest)

# need CGAL 4.7
# cmake -DCGAL_DIR=~/workspace/CGAL-4.7-Ic-41 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )
 
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      include_directories (BEFORE "../../include")

      add_executable ( AlphaComplexUT Alpha_complex_unit_test.cpp )
      target_link_libraries(AlphaComplexUT  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

      add_executable ( cerr cerr.cpp )

      # Do not forget to copy test files in current binary dir
      file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

      add_test(AlphaComplexUT ${CMAKE_CURRENT_BINARY_DIR}/AlphaComplexUT
      # XML format for Jenkins xUnit plugin 
      --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/AlphaComplexUT.xml --log_level=test_suite --report_level=no)

    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha complex feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha complex feature. Version 4.6.0 is required.")
  endif ()
endif()

