cmake_minimum_required(VERSION 2.6)
project(Hasse_diagram)

add_executable ( Hasse_diagram_example Hasse_diagram_example.cpp )
if (TBB_FOUND)
  target_link_libraries(Hasse_diagram_example ${TBB_LIBRARIES})
endif()

add_executable ( Hasse_diagram_from_cubical_complex Hasse_diagram_from_cubical_complex.cpp )
if (TBB_FOUND)
  target_link_libraries(Hasse_diagram_from_cubical_complex ${TBB_LIBRARIES})
endif()

add_executable ( Hasse_diagram_from_simplex_tree Hasse_diagram_from_simplex_tree.cpp )
if (TBB_FOUND)
  target_link_libraries(Hasse_diagram_from_simplex_tree ${TBB_LIBRARIES})
endif()

add_executable ( Insertion_and_removal_operations_on_Hasse_diagrams Insertion_and_removal_operations_on_Hasse_diagrams.cpp )
if (TBB_FOUND)
  target_link_libraries(Insertion_and_removal_operations_on_Hasse_diagrams ${TBB_LIBRARIES})
endif()

install(TARGETS Hasse_diagram_example DESTINATION bin)
