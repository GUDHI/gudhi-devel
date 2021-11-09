project(Discrete_morse_theory_examples)

add_executable(simple_discrete_morse_theory simple_discrete_morse_theory.cpp)
if (TBB_FOUND)
  target_link_libraries(simple_discrete_morse_theory ${TBB_LIBRARIES})
endif()