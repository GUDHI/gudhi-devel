add_executable(ripser ripser.cc)
install(TARGETS ripser DESTINATION bin)
add_test(NAME ripser_lower COMMAND $<TARGET_FILE:ripser> --format lower --dim 2 --threshold .7 "${CMAKE_SOURCE_DIR}/data/distance_matrix/tore3D_1307_distance_matrix.csv")
