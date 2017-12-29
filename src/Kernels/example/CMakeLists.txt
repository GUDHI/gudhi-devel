cmake_minimum_required(VERSION 2.6)
project(Kernels_examples)

add_executable ( BasicEx kernel_basic_example.cpp )

if (TBB_FOUND)
  target_link_libraries(BasicEx ${TBB_LIBRARIES})
endif()

add_test(NAME Kernels_example_basicex COMMAND $<TARGET_FILE:BasicEx>
    "")
