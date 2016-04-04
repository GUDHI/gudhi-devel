# CMake generated Testfile for 
# Source directory: /home/frg/Bureau/mWorkingDirectory/src/Persistent_cohomology/test
# Build directory: /home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(PersistentCohomologyUT "/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/test/PersistentCohomologyUT" "/home/frg/Bureau/mWorkingDirectory/src/Persistent_cohomology/test/simplex_tree_file_for_unit_test.txt" "--log_format=XML" "--log_sink=/home/frg/Bureau/mWorkingDirectory/PersistentCohomologyUT.xml" "--log_level=test_suite" "--report_level=no")
add_test(PersistentCohomologyMultiFieldUT "/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/test/PersistentCohomologyMultiFieldUT" "/home/frg/Bureau/mWorkingDirectory/src/Persistent_cohomology/test/simplex_tree_file_for_multi_field_unit_test.txt" "--log_format=XML" "--log_sink=/home/frg/Bureau/mWorkingDirectory/PersistentCohomologyMultiFieldUT.xml" "--log_level=test_suite" "--report_level=no")
