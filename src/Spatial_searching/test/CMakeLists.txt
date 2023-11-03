project(Spatial_searching_tests)

if(CGAL_FOUND AND TARGET Eigen3::Eigen)
  include(GUDHI_boost_test)
  add_executable( Spatial_searching_test_Kd_tree_search test_Kd_tree_search.cpp )
  target_link_libraries(Spatial_searching_test_Kd_tree_search CGAL::CGAL Eigen3::Eigen)
  gudhi_add_boost_test(Spatial_searching_test_Kd_tree_search)
endif ()
