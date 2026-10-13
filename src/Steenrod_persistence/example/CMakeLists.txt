add_executable_with_targets(Steenrod_persistence_example_rp2 rp2_steenrod.cpp TBB::tbb)

add_test(NAME Steenrod_persistence_example_rp2
         COMMAND $<TARGET_FILE:Steenrod_persistence_example_rp2>)
