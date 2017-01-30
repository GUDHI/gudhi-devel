cmake_minimum_required(VERSION 2.6)
project(Nerve_GIC_examples)

add_executable ( GIC simple_GIC.cpp )
target_link_libraries(GIC ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(GIC ${TBB_LIBRARIES})
endif()