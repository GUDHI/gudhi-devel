cmake_minimum_required(VERSION 2.6)
project(Hasse_diagram)

add_executable ( Hasse_diagram_from_file Hasse_diagram_from_file.cpp )
if (TBB_FOUND)
  target_link_libraries(Hasse_diagram_from_file ${TBB_LIBRARIES})
endif()

install(TARGETS Hasse_diagram_from_file DESTINATION bin)
