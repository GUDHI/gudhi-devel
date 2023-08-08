project(Zigzag_benchmark)

find_package(benchmark REQUIRED)

add_executable(Oscillating_rips_benchmark Oscillating_rips_benchmark.cpp)
target_link_libraries(Oscillating_rips_benchmark benchmark::benchmark)
if(TARGET TBB::tbb)
  target_link_libraries(Oscillating_rips_benchmark TBB::tbb)
endif()

add_executable(Comp_benchmark Oscillating_rips_comp_benchmark.cpp ./ext_zz/fzz/fzz.cpp)
target_link_libraries(Comp_benchmark benchmark::benchmark)
target_include_directories(Comp_benchmark PUBLIC . ./ext_zz)
if(TARGET TBB::tbb)
  target_link_libraries(Comp_benchmark TBB::tbb)
endif()
target_compile_options(Comp_benchmark PUBLIC "-fopenmp")
target_link_options(Comp_benchmark PUBLIC "-fopenmp")

