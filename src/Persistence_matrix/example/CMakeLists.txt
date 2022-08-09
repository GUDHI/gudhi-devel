project(Matrix_examples)

set(CMAKE_CXX_STANDARD 17)

add_executable(Matrix_examples_test comp_test.cpp)
add_test(NAME Matrix_examples_test COMMAND $<TARGET_FILE:Matrix_examples_test>)
