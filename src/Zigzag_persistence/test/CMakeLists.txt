project(Zigzag_persistence_tests)

include(GUDHI_boost_test)

add_executable ( Zigzag_persistence_unit_test zigzag_persistence_unit_test.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(Zigzag_persistence_unit_test TBB::tbb)
endif()

gudhi_add_boost_test(Zigzag_persistence_unit_test)

