cmake_minimum_required(VERSION 2.6)
project(GUDHIExPersCohom)

# NEED TO FIND BOOST NEEDED COMPONENTS TO LINK WITH
find_package(Boost 1.45.0 COMPONENTS system program_options)
message("Boost_lib = ${Boost_LIBRARIES}")

add_executable ( rips_persistence rips_persistence.cpp )
target_link_libraries(rips_persistence ${Boost_LIBRARIES})
add_executable ( rips_persistence_from_alpha_shape_3 rips_persistence_from_alpha_shape_3.cpp )
target_link_libraries(rips_persistence_from_alpha_shape_3 ${Boost_LIBRARIES})

if(GMPXX_FOUND AND GMP_FOUND)
   message("GMPXX_LIBRARIES = ${GMPXX_LIBRARIES}")
   message("GMP_LIBRARIES = ${GMP_LIBRARIES}")

   add_executable ( rips_multifield_persistence rips_multifield_persistence.cpp )
   target_link_libraries(rips_multifield_persistence  ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

   add_executable ( performance_rips_persistence performance_rips_persistence.cpp )
   target_link_libraries(performance_rips_persistence  ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
endif()
