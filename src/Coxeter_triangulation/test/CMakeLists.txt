cmake_minimum_required(VERSION 2.6)
project(Coxeter_triangulation_test)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
  add_executable ( Coxeter_triangulation_coxeter_test coxeter_test.cpp )
  if (TBB_FOUND)
    target_link_libraries(Coxeter_triangulation_coxeter_test ${TBB_LIBRARIES})
  endif()
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
