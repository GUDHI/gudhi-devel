project(Alpha_complex_tests)

include(GUDHI_boost_test)
if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)

  # Do not forget to copy test files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  add_executable ( Alpha_complex_test_unit Alpha_complex_unit_test.cpp )
  target_link_libraries(Alpha_complex_test_unit ${CGAL_LIBRARY})
  add_executable ( Delaunay_complex_test_unit Delaunay_complex_unit_test.cpp )
  target_link_libraries(Delaunay_complex_test_unit ${CGAL_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Alpha_complex_test_unit ${TBB_LIBRARIES})
    target_link_libraries(Delaunay_complex_test_unit ${TBB_LIBRARIES})
  endif()

  gudhi_add_boost_test(Alpha_complex_test_unit)
  gudhi_add_boost_test(Delaunay_complex_test_unit)

endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)

if (NOT CGAL_VERSION VERSION_LESS 4.11.0)

  add_executable ( Alpha_complex_3d_test_unit Alpha_complex_3d_unit_test.cpp )
  target_link_libraries(Alpha_complex_3d_test_unit ${CGAL_LIBRARY})
  add_executable ( Weighted_alpha_complex_3d_test_unit Weighted_alpha_complex_3d_unit_test.cpp )
  target_link_libraries(Weighted_alpha_complex_3d_test_unit ${CGAL_LIBRARY})
  add_executable ( Periodic_alpha_complex_3d_test_unit Periodic_alpha_complex_3d_unit_test.cpp )
  target_link_libraries(Periodic_alpha_complex_3d_test_unit ${CGAL_LIBRARY})
  add_executable ( Weighted_periodic_alpha_complex_3d_test_unit Weighted_periodic_alpha_complex_3d_unit_test.cpp )
  target_link_libraries(Weighted_periodic_alpha_complex_3d_test_unit ${CGAL_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Alpha_complex_3d_test_unit ${TBB_LIBRARIES})
    target_link_libraries(Weighted_alpha_complex_3d_test_unit ${TBB_LIBRARIES})
    target_link_libraries(Periodic_alpha_complex_3d_test_unit ${TBB_LIBRARIES})
    target_link_libraries(Weighted_periodic_alpha_complex_3d_test_unit ${TBB_LIBRARIES})
  endif()

  gudhi_add_boost_test(Alpha_complex_3d_test_unit)
  gudhi_add_boost_test(Weighted_alpha_complex_3d_test_unit)
  gudhi_add_boost_test(Periodic_alpha_complex_3d_test_unit)
  gudhi_add_boost_test(Weighted_periodic_alpha_complex_3d_test_unit)

endif (NOT CGAL_VERSION VERSION_LESS 4.11.0)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 5.1.0)
  add_executable ( Alpha_kernel_d_test_unit Alpha_kernel_d_unit_test.cpp )
  target_link_libraries(Alpha_kernel_d_test_unit ${CGAL_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Alpha_kernel_d_test_unit ${TBB_LIBRARIES})
  endif()
  gudhi_add_boost_test(Alpha_kernel_d_test_unit)
  
  add_executable ( Weighted_alpha_complex_test_unit Weighted_alpha_complex_unit_test.cpp )
  target_link_libraries(Weighted_alpha_complex_test_unit ${CGAL_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Weighted_alpha_complex_test_unit ${TBB_LIBRARIES})
  endif()
  gudhi_add_boost_test(Weighted_alpha_complex_test_unit)

endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 5.1.0)