project(Zigzag_persistence_examples)

add_executable ( Zigzag_persistence_example_simple_zigzag_filtration example_simple_zigzag_filtration.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(Zigzag_persistence_example_simple_zigzag_filtration TBB::tbb)
endif()
add_test(NAME Zigzag_persistence_example_simple_zigzag_filtration COMMAND $<TARGET_FILE:Zigzag_persistence_example_simple_zigzag_filtration>)

add_executable ( Zigzag_persistence_example_zzfiltration_from_file example_zzfiltration_from_file.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(Zigzag_persistence_example_zzfiltration_from_file TBB::tbb)
endif()
file(COPY "zigzag_filtration_example.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
add_test(NAME Zigzag_persistence_example_zzfiltration_from_file COMMAND $<TARGET_FILE:Zigzag_persistence_example_zzfiltration_from_file> "${CMAKE_CURRENT_BINARY_DIR}/zigzag_filtration_example.txt")

add_executable ( Zigzag_persistence_example_oscillating_rips_persistence example_oscillating_rips_persistence.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(Zigzag_persistence_example_oscillating_rips_persistence TBB::tbb)
endif()
add_test(NAME Zigzag_persistence_example_oscillating_rips_persistence COMMAND $<TARGET_FILE:Zigzag_persistence_example_oscillating_rips_persistence> "2" "3" "5" "10")

add_executable ( comp comparison_for_tests.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(comp TBB::tbb)
endif()




