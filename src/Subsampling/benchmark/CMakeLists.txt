project(Subsampling_benchmark)

add_executable ( choose_n_farthest_points choose_n_farthest_points.cpp )
add_test(NAME Compare_choose_n_farthest_points COMMAND $<TARGET_FILE:choose_n_farthest_points>)
