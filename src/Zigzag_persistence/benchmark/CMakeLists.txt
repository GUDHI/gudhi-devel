project(Zigzag_benchmark)

find_package(benchmark REQUIRED)

add_executable(Zigzag_benchmark Zigzag_benchmark.cpp)
target_link_libraries(Zigzag_benchmark benchmark::benchmark)
target_include_directories(Zigzag_benchmark PUBLIC . ../example/ext_zz)

add_executable(Zigzag_old_benchmark Zigzag_old_benchmark.cpp)
target_link_libraries(Zigzag_old_benchmark benchmark::benchmark)
target_include_directories(Zigzag_old_benchmark PUBLIC . ../example/ext_zz)

add_executable(Zigzag_other_benchmark Zigzag_other_benchmark.cpp ../example/ext_zz/fzz/fzz.cpp)
target_link_libraries(Zigzag_other_benchmark benchmark::benchmark)
target_include_directories(Zigzag_other_benchmark PUBLIC . ../example/ext_zz)
target_compile_options(Zigzag_other_benchmark PUBLIC "-fopenmp")
target_link_options(Zigzag_other_benchmark PUBLIC "-fopenmp")
