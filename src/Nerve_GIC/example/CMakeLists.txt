cmake_minimum_required(VERSION 2.6)
project(Nerve_GIC_examples)

add_executable ( Nerve Nerve.cpp )
target_link_libraries(Nerve  ${Boost_SYSTEM_LIBRARY})

add_executable ( CoordGIC CoordGIC.cpp )
target_link_libraries(CoordGIC  ${Boost_SYSTEM_LIBRARY})

add_executable ( FuncGIC FuncGIC.cpp )
target_link_libraries(FuncGIC  ${Boost_SYSTEM_LIBRARY})

add_executable ( VoronoiGIC VoronoiGIC.cpp )
target_link_libraries(VoronoiGIC  ${Boost_SYSTEM_LIBRARY})

#[[file(COPY KeplerMapperVisuFromTxtFile.py km.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)]]

if (TBB_FOUND)
  target_link_libraries(Nerve ${TBB_LIBRARIES})
  target_link_libraries(CoordGIC ${TBB_LIBRARIES})
  target_link_libraries(FuncGIC ${TBB_LIBRARIES})
  target_link_libraries(VoronoiGIC ${TBB_LIBRARIES})
endif()