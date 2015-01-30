cmake_minimum_required(VERSION 2.6)
project(GUDHIskbl)


add_executable(SkeletonBlockerFromSimplices Skeleton_blocker_from_simplices.cpp)

add_executable(SkeletonBlockerIteration Skeleton_blocker_iteration.cpp)
target_link_libraries(SkeletonBlockerIteration ${Boost_TIMER_LIBRARY} ${Boost_SYSTEM_LIBRARY})	
