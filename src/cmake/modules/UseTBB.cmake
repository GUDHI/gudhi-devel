# This module setups the compiler for using TBB library.
# It assumes that find_package(TBB) was already called.

include_directories ( ${TBB_INCLUDE_DIRS} )
link_directories( ${TBB_LIBRARY_DIRS} )
add_definitions( -DNOMINMAX -DGUDHI_USE_TBB )
