add_executable ( Coxeter_triangulation_implicit_manifold_example implicit_manifold.cpp )
add_executable ( Coxeter_triangulation_lowdim_manifold_example coxeter_lowdim_manif_time.cpp )
add_executable ( Coxeter_triangulation_chair_cox_fkt_example coxeter_chair_cox_fkt.cpp )
target_link_libraries( Coxeter_triangulation_implicit_manifold_example ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )

if (TBB_FOUND)
  target_link_libraries( Coxeter_triangulation_implicit_manifold_example ${TBB_LIBRARIES} )
endif()
