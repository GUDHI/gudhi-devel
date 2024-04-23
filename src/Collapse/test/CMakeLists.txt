add_executable_with_targets(Collapse_test_unit collapse_unit_test.cpp Eigen3::Eigen TBB::tbb)

if (TARGET Eigen3::Eigen)
  include(GUDHI_boost_test)
  gudhi_add_boost_test(Collapse_test_unit)
endif()