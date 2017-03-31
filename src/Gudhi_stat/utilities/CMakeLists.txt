cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT)



#some stat libraries:
add_executable ( permutation_test permutation_test.cpp )
target_link_libraries(permutation_test ${Boost_SYSTEM_LIBRARY})

#add_executable ( topological_process topological_process.cpp )
#target_link_libraries(topological_process ${Boost_SYSTEM_LIBRARY})

#add_executable ( topological_process_2 topological_process_2.cpp )
#target_link_libraries(topological_process_2 ${Boost_SYSTEM_LIBRARY})

add_executable ( Hausdorff_subsampling Hausdorff_subsampling.cpp )
if (TBB_FOUND)
target_link_libraries(Hausdorff_subsampling ${TBB_LIBRARIES})
endif(TBB_FOUND)
target_link_libraries(Hausdorff_subsampling ${Boost_SYSTEM_LIBRARY})


add_executable ( Landscape_bootstrap Landscape_bootstrap.cpp )
if (TBB_FOUND)
target_link_libraries(Landscape_bootstrap ${TBB_LIBRARIES})
endif(TBB_FOUND)
target_link_libraries(Landscape_bootstrap ${Boost_SYSTEM_LIBRARY})


add_executable ( Multiplicative_bootstrap Multiplicative_bootstrap.cpp )
if (TBB_FOUND)
target_link_libraries(Multiplicative_bootstrap ${TBB_LIBRARIES})
endif(TBB_FOUND)
target_link_libraries(Multiplicative_bootstrap ${Boost_SYSTEM_LIBRARY})




#just for Pawel's tests:
#add_executable ( NN_classification NN_classification.cpp )
#target_link_libraries(NN_classification ${Boost_SYSTEM_LIBRARY})
#
#add_executable ( compute_distance_between_vectors compute_distance_between_vectors.cpp )
#target_link_libraries(compute_distance_between_vectors ${Boost_SYSTEM_LIBRARY})


if(CGAL_FOUND)
  add_executable ( sliding_window_embedding sliding_window_embedding.cpp )
  target_link_libraries(sliding_window_embedding ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(sliding_window_embedding ${TBB_LIBRARIES})
  endif(TBB_FOUND)
endif(CGAL_FOUND)
