# A function to add a new module in GUDHI

set(GUDHI_MODULES "")
function(add_gudhi_module file_path)
    set(GUDHI_MODULES ${GUDHI_MODULES} ${file_path} PARENT_SCOPE)
endfunction(add_gudhi_module)

# Add your new module in the list, order is not important

add_gudhi_module(common)
add_gudhi_module(Alpha_complex)
add_gudhi_module(Bitmap_cubical_complex)
add_gudhi_module(Bottleneck_distance)
add_gudhi_module(Contraction)
add_gudhi_module(Hasse_complex)
add_gudhi_module(Persistent_cohomology)
add_gudhi_module(Rips_complex)
add_gudhi_module(Simplex_tree)
add_gudhi_module(Skeleton_blocker)
add_gudhi_module(Spatial_searching)
add_gudhi_module(Subsampling)
add_gudhi_module(Tangential_complex)
add_gudhi_module(Witness_complex)

# message("++ GUDHI_MODULES list is:\"${GUDHI_MODULES}\"")

if (WITH_GUDHI_BENCHMARK)
  set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};benchmark")
endif()
if (WITH_GUDHI_EXAMPLE)
  set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};example")
endif()
if (NOT WITHOUT_GUDHI_TEST)
    set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};test")
endif()
if (NOT WITHOUT_GUDHI_UTILITIES)
  set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};utilities")
endif()

message("++ GUDHI_SUB_DIRECTORIES list is:\"${GUDHI_SUB_DIRECTORIES}\"")



