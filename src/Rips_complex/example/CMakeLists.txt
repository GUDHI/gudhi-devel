cmake_minimum_required(VERSION 2.6)
project(Rips_complex_examples)

# Point cloud
add_executable ( ripsoffreader example_rips_complex_from_off_file.cpp )
target_link_libraries(ripsoffreader  ${Boost_SYSTEM_LIBRARY})

add_executable ( oneskeletonripspoints example_one_skeleton_rips_from_points.cpp )
target_link_libraries(oneskeletonripspoints  ${Boost_SYSTEM_LIBRARY})

# Distance matrix
add_executable ( oneskeletonripsdistance example_one_skeleton_rips_from_distance_matrix.cpp )
target_link_libraries(oneskeletonripsdistance  ${Boost_SYSTEM_LIBRARY})

add_executable ( ripscsvdistancereader example_rips_complex_from_csv_distance_matrix_file.cpp )
target_link_libraries(ripscsvdistancereader  ${Boost_SYSTEM_LIBRARY})

if (TBB_FOUND)
  target_link_libraries(ripsoffreader ${TBB_LIBRARIES})
  target_link_libraries(oneskeletonripspoints ${TBB_LIBRARIES})
  target_link_libraries(oneskeletonripsdistance ${TBB_LIBRARIES})
  target_link_libraries(ripscsvdistancereader ${TBB_LIBRARIES})
endif()

add_test(oneskeletonripspoints ${CMAKE_CURRENT_BINARY_DIR}/oneskeletonripspoints 12.0)

add_test(oneskeletonripsdistance ${CMAKE_CURRENT_BINARY_DIR}/oneskeletonripspoints 1.0)

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
add_test(ripsoffreader_doc_12_1 ${CMAKE_CURRENT_BINARY_DIR}/ripsoffreader alphacomplexdoc.off 12.0 1 ${CMAKE_CURRENT_BINARY_DIR}/ripsoffreader_result_12_1.txt)
add_test(ripsoffreader_doc_12_3 ${CMAKE_CURRENT_BINARY_DIR}/ripsoffreader alphacomplexdoc.off 12.0 3 ${CMAKE_CURRENT_BINARY_DIR}/ripsoffreader_result_12_3.txt)

file(COPY "${CMAKE_SOURCE_DIR}/data/distance_matrix/full_square_distance_matrix.csv" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
add_test(ripscsvdistancereader_doc_1_1 ${CMAKE_CURRENT_BINARY_DIR}/ripscsvdistancereader full_square_distance_matrix.csv 1.0 1 ${CMAKE_CURRENT_BINARY_DIR}/ripscsvreader_result_1_1.txt)
add_test(ripscsvdistancereader_doc_1_3 ${CMAKE_CURRENT_BINARY_DIR}/ripscsvdistancereader full_square_distance_matrix.csv 1.0 3 ${CMAKE_CURRENT_BINARY_DIR}/ripscsvreader_result_1_3.txt)


if (DIFF_PATH)
  # Do not forget to copy test results files in current binary dir
  file(COPY "one_skeleton_rips_points_for_doc.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY "full_skeleton_rips_points_for_doc.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY "one_skeleton_rips_distance_for_doc.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY "full_skeleton_rips_distance_for_doc.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  add_test(ripsoffreader_doc_12_1_diff_files ${DIFF_PATH} ${CMAKE_CURRENT_BINARY_DIR}/ripsoffreader_result_12_1.txt ${CMAKE_CURRENT_BINARY_DIR}/one_skeleton_rips_points_for_doc.txt)
  add_test(ripsoffreader_doc_12_3_diff_files ${DIFF_PATH} ${CMAKE_CURRENT_BINARY_DIR}/ripsoffreader_result_12_3.txt ${CMAKE_CURRENT_BINARY_DIR}/full_skeleton_rips_points_for_doc.txt)
  add_test(ripscsvreader_doc_1_1_diff_files ${DIFF_PATH} ${CMAKE_CURRENT_BINARY_DIR}/ripscsvreader_result_1_1.txt ${CMAKE_CURRENT_BINARY_DIR}/one_skeleton_rips_distance_for_doc.txt)
  add_test(ripscsvreader_doc_1_3_diff_files ${DIFF_PATH} ${CMAKE_CURRENT_BINARY_DIR}/ripscsvreader_result_1_3.txt ${CMAKE_CURRENT_BINARY_DIR}/full_skeleton_rips_distance_for_doc.txt)
endif()
