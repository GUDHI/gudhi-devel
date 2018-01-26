cmake_minimum_required(VERSION 2.6)
project(Nerve_GIC_examples)

if (NOT CGAL_VERSION VERSION_LESS 4.8.1)

  add_executable ( Nerve Nerve.cpp )
  add_executable ( VoronoiGIC VoronoiGIC.cpp )

  if (TBB_FOUND)
    target_link_libraries(Nerve ${TBB_LIBRARIES})
    target_link_libraries(VoronoiGIC ${TBB_LIBRARIES})
  endif()

  file(COPY KeplerMapperVisuFromTxtFile.py km.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  add_test(NAME Nerve_GIC_utilities_nerve COMMAND $<TARGET_FILE:Nerve>
      "${CMAKE_SOURCE_DIR}/data/points/human.off" "2" "10" "0.3")

  add_test(NAME Nerve_GIC_utilities_VoronoiGIC COMMAND $<TARGET_FILE:VoronoiGIC>
      "${CMAKE_SOURCE_DIR}/data/points/human.off" "100")

endif (NOT CGAL_VERSION VERSION_LESS 4.8.1)
