project(Zigzag_persistence_examples)

add_executable ( Zigzag_persistence_example_simple_zigzag_filtration example_simple_zigzag_filtration.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(Zigzag_persistence_example_simple_zigzag_filtration TBB::tbb)
endif()
add_test(NAME Zigzag_persistence_example_simple_zigzag_filtration COMMAND $<TARGET_FILE:Zigzag_persistence_example_simple_zigzag_filtration>)

add_executable ( comp comparison_for_tests.cpp ./ext_zz/fzz/fzz.cpp )
target_include_directories(comp PUBLIC ./ext_zz)
target_compile_options(comp PUBLIC "-fopenmp")
target_link_options(comp PUBLIC "-fopenmp")

add_executable ( Zigzag_persistence_example_rips_zigzag_filtration example_rips_zigzag_filtration.cpp )
target_include_directories(Zigzag_persistence_example_rips_zigzag_filtration PUBLIC ./ext_zz)
target_compile_options(Zigzag_persistence_example_rips_zigzag_filtration PUBLIC "-fopenmp")
target_link_options(Zigzag_persistence_example_rips_zigzag_filtration PUBLIC "-fopenmp")

# add_executable ( rips example_rips_zigzag_filtration.cpp )
# target_include_directories(rips PUBLIC ./ext_zz)
# target_compile_options(rips PUBLIC "-fopenmp")
# target_link_options(rips PUBLIC "-fopenmp")

# add_executable ( rips_old example_rips_zigzag_filtration.cpp )
# target_include_directories(rips_old PUBLIC ./ext_zz)
# target_compile_options(rips_old PUBLIC "-fopenmp")
# target_link_options(rips_old PUBLIC "-fopenmp")


