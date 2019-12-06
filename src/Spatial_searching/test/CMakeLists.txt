project(Spatial_searching_tests)

if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  include(GUDHI_boost_test)
  add_executable( Spatial_searching_test_Kd_tree_search test_Kd_tree_search.cpp )
  target_link_libraries(Spatial_searching_test_Kd_tree_search ${CGAL_LIBRARY})
  gudhi_add_boost_test(Spatial_searching_test_Kd_tree_search)
endif ()
