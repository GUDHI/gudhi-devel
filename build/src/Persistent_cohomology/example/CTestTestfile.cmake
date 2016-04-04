# CMake generated Testfile for 
# Source directory: /home/frg/Bureau/mWorkingDirectory/src/Persistent_cohomology/example
# Build directory: /home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/example
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(rips_persistence_2 "/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/example/rips_persistence" "/home/frg/Bureau/mWorkingDirectory/data/points/Kl.txt" "-r" "0.25" "-d" "3" "-p" "2" "-m" "100")
add_test(rips_persistence_3 "/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/example/rips_persistence" "/home/frg/Bureau/mWorkingDirectory/data/points/Kl.txt" "-r" "0.25" "-d" "3" "-p" "3" "-m" "100")
add_test(persistence_from_file_3_2_0 "/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/example/persistence_from_file" "/home/frg/Bureau/mWorkingDirectory/data/points/bunny_5000.st" "-p" "2" "-m" "0")
add_test(persistence_from_file_3_3_100 "/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/example/persistence_from_file" "/home/frg/Bureau/mWorkingDirectory/data/points/bunny_5000.st" "-p" "3" "-m" "100")
add_test(rips_multifield_persistence_2_3 "/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/example/rips_multifield_persistence" "/home/frg/Bureau/mWorkingDirectory/data/points/Kl.txt" "-r" "0.25" "-d" "3" "-p" "2" "-q" "3" "-m" "100")
add_test(rips_multifield_persistence_2_71 "/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/example/rips_multifield_persistence" "/home/frg/Bureau/mWorkingDirectory/data/points/Kl.txt" "-r" "0.25" "-d" "3" "-p" "2" "-q" "71" "-m" "100")
add_test(alpha_shapes_persistence_2_0_5 "/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/example/alpha_shapes_persistence" "/home/frg/Bureau/mWorkingDirectory/data/points/bunny_5000" "2" "0.5")
