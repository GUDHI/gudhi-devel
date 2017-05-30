# A function to add a new module in GUDHI

set(GUDHI_MODULES "")
function(add_gudhi_module file_path)
    set(GUDHI_MODULES ${GUDHI_MODULES} ${file_path} PARENT_SCOPE)
endfunction(add_gudhi_module)

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



