cmake_minimum_required(VERSION 2.6)
project(Nerve_GIC_examples)

add_executable ( Nerve Nerve.cpp )
target_link_libraries(Nerve  ${Boost_SYSTEM_LIBRARY})

add_executable ( GIC GIC.cpp )
target_link_libraries(GIC  ${Boost_SYSTEM_LIBRARY})

add_executable ( MapperDeltaCoord MapperDeltaCoord.cpp )
target_link_libraries(MapperDeltaCoord  ${Boost_SYSTEM_LIBRARY})

add_executable ( MapperDeltaFunc MapperDeltaFunc.cpp )
target_link_libraries(MapperDeltaFunc  ${Boost_SYSTEM_LIBRARY})

add_executable ( GICvoronoi GICvoronoi.cpp )
target_link_libraries(GICvoronoi  ${Boost_SYSTEM_LIBRARY})

file(COPY visu.py km.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

if (TBB_FOUND)
  target_link_libraries(Nerve ${TBB_LIBRARIES})
  target_link_libraries(GIC ${TBB_LIBRARIES})
  target_link_libraries(MapperDeltaCoord ${TBB_LIBRARIES})
  target_link_libraries(MapperDeltaFunc ${TBB_LIBRARIES})
  target_link_libraries(GICvoronoi ${TBB_LIBRARIES})
endif()