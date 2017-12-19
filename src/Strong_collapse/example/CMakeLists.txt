cmake_minimum_required(VERSION 2.6)
project(Strong_collapse)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
  add_executable (example_usage example_usage.cpp )
  target_link_libraries(example_usage ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(example_usage ${TBB_LIBRARIES})
  endif()

  add_executable (example_usage2 example_usage2.cpp )
  target_link_libraries(example_usage2 ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(example_usage2 ${TBB_LIBRARIES})
  endif()

  add_executable (alpha1 alpha1.cpp )
  target_link_libraries(alpha1 ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(alpha1 ${TBB_LIBRARIES})
  endif()

  add_executable (alpha2 alpha2.cpp )
  target_link_libraries(alpha2 ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(alpha2 ${TBB_LIBRARIES})
  endif()

  add_executable (rips1 rips1.cpp )
  target_link_libraries(rips1 ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(rips1 ${TBB_LIBRARIES})
  endif()

  add_executable (rips2 rips2.cpp )
  target_link_libraries(rips2 ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(rips2 ${TBB_LIBRARIES})
  endif()
  
  add_executable (sparseCollapseRips sparseCollapseRips.cpp )
  target_link_libraries(sparseCollapseRips ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(sparseCollapseRips ${TBB_LIBRARIES})
  endif()
  
  add_executable (manualSimplicialComplexColl manualSimplicialComplexColl.cpp )
  target_link_libraries(manualSimplicialComplexColl ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(manualSimplicialComplexColl ${TBB_LIBRARIES})
  endif()

  add_executable (alphaSparseCollapse alphaSparseCollapse.cpp )
  target_link_libraries(alphaSparseCollapse ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
   target_link_libraries(alphaSparseCollapse ${TBB_LIBRARIES})
  endif()	

	
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)

