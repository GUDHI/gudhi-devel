project(Toplex_map_tests)

include(GUDHI_boost_test)

add_executable( Toplex_map_unit_test toplex_map_unit_test.cpp )
gudhi_add_boost_test(Toplex_map_unit_test)

add_executable( Lazy_toplex_map_unit_test lazy_toplex_map_unit_test.cpp )
gudhi_add_boost_test(Lazy_toplex_map_unit_test)
