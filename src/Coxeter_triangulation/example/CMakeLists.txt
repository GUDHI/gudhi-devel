add_executable ( Coxeter_triangulation_implicit_manifold_example implicit_manifold.cpp )
target_link_libraries( Coxeter_triangulation_implicit_manifold_example ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )

if (TBB_FOUND)
  target_link_libraries( Coxeter_triangulation_implicit_manifold_example ${TBB_LIBRARIES} )
endif()
