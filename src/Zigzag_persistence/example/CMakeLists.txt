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

# add_executable ( comp comparison_for_tests.cpp ./ext_zz/fzz/fzz.cpp )
# target_include_directories(comp PUBLIC ./ext_zz)
# target_compile_options(comp PUBLIC "-fopenmp")
# target_link_options(comp PUBLIC "-fopenmp")

# add_executable ( Zigzag_persistence_example_rips_zigzag_filtration example_rips_zigzag_filtration.cpp )
# target_include_directories(Zigzag_persistence_example_rips_zigzag_filtration PUBLIC ./ext_zz)
# target_compile_options(Zigzag_persistence_example_rips_zigzag_filtration PUBLIC "-fopenmp")
# target_link_options(Zigzag_persistence_example_rips_zigzag_filtration PUBLIC "-fopenmp")

# add_executable ( rips example_rips_zigzag_filtration.cpp )
# target_include_directories(rips PUBLIC ./ext_zz)
# target_compile_options(rips PUBLIC "-fopenmp")
# target_link_options(rips PUBLIC "-fopenmp")

# add_executable ( rips_old example_rips_zigzag_filtration.cpp )
# target_include_directories(rips_old PUBLIC ./ext_zz)
# target_compile_options(rips_old PUBLIC "-fopenmp")
# target_link_options(rips_old PUBLIC "-fopenmp")


