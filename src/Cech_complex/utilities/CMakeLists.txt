cmake_minimum_required(VERSION 2.6)
project(Cech_complex_utilities)

add_executable(cech_persistence cech_persistence.cpp)
target_link_libraries(cech_persistence ${Boost_PROGRAM_OPTIONS_LIBRARY})

if (TBB_FOUND)
  target_link_libraries(cech_persistence ${TBB_LIBRARIES})
endif()

add_test(NAME Cech_complex_utility_from_rips_on_tore_3D COMMAND $<TARGET_FILE:cech_persistence>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3")

install(TARGETS cech_persistence DESTINATION bin)
