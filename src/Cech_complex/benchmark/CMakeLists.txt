cmake_minimum_required(VERSION 2.6)
project(Cech_complex_benchmark)

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_executable(cech_complex_benchmark cech_complex_benchmark.cpp)
target_link_libraries(cech_complex_benchmark ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

if (TBB_FOUND)
  target_link_libraries(cech_complex_benchmark ${TBB_LIBRARIES})
endif()
