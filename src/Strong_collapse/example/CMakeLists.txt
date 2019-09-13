project(Strong_collapse_examples)

add_executable(strong_collapse_from_points strong_collapse_from_points.cpp)
if (TBB_FOUND)
  target_link_libraries(strong_collapse_from_points ${TBB_LIBRARIES})
endif()

add_executable(strong_collapse_rips_persistence_step_by_step strong_collapse_rips_persistence_step_by_step.cpp)
target_link_libraries(strong_collapse_rips_persistence_step_by_step ${Boost_PROGRAM_OPTIONS_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(strong_collapse_rips_persistence_step_by_step ${TBB_LIBRARIES})
endif()
