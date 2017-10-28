cmake_minimum_required(VERSION 2.6)
project(Persistent_cohomology_benchmark)

  
if(GMP_FOUND)
   if(GMPXX_FOUND)
      add_executable ( performance_rips_persistence EXCLUDE_FROM_ALL performance_rips_persistence.cpp )
      target_link_libraries(performance_rips_persistence ${Boost_PROGRAM_OPTIONS_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
      if (TBB_FOUND)
        target_link_libraries(performance_rips_persistence ${TBB_LIBRARIES})
      endif(TBB_FOUND)
      file(COPY "${CMAKE_SOURCE_DIR}/data/points/Kl.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
   endif(GMPXX_FOUND)
endif(GMP_FOUND)
