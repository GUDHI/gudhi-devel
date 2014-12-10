cmake_minimum_required(VERSION 2.6)
project(GUDHITestSimplexTree)

# NEED TO FIND BOOST NEEDED COMPONENTS TO LINK WITH
find_package(Boost 1.45.0 COMPONENTS system unit_test_framework)
message("Boost_lib = ${Boost_LIBRARIES}")

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} --coverage")
endif()

include_directories(${Boost_INCLUDE_DIRS})
add_executable ( persistent_cohomology_unit_test persistent_cohomology_unit_test.cpp )
target_link_libraries(persistent_cohomology_unit_test ${Boost_LIBRARIES})

# Unitary tests
add_test(persistent_cohomology_unit_test ${CMAKE_CURRENT_BINARY_DIR}/persistent_cohomology_unit_test)

if(GMPXX_FOUND AND GMP_FOUND)
   add_executable ( persistent_cohomology_unit_test_multi_field persistent_cohomology_unit_test_multi_field.cpp )
   target_link_libraries(persistent_cohomology_unit_test_multi_field  ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

   # Unitary tests
   add_test(persistent_cohomology_unit_test_multi_field ${CMAKE_CURRENT_BINARY_DIR}/persistent_cohomology_unit_test_multi_field)
endif()

if (LCOV_PATH)
  # Lcov code coverage of unitary test
  add_test(persistent_cohomology_unit_test_coverage.info ${CMAKE_SOURCE_DIR}/scripts/check_code_coverage.sh ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology)
endif()

if (PYTHON_PATH)
  # Cpplint tests on coding style
  add_test(Persistent_cohomology.h.cpplint ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/include/gudhi/Persistent_cohomology.h)
  add_test(Field_Zp.h.cpplint ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/include/gudhi/Persistent_cohomology/Field_Zp.h)
  add_test(Multi_field.h.cpplint ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/include/gudhi/Persistent_cohomology/Multi_field.h)
  add_test(Persistent_cohomology_column.h.cpplint ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/include/gudhi/Persistent_cohomology/Persistent_cohomology_column.h)
endif()