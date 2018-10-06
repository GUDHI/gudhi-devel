# A function to add a new module in GUDHI

set(GUDHI_MODULES_FULL_LIST "")
function(add_gudhi_module file_path)
  option("WITH_MODULE_GUDHI_${file_path}" "Activate/desactivate ${file_path} compilation and installation" ON)
  if (WITH_MODULE_GUDHI_${file_path})
    set(GUDHI_MODULES ${GUDHI_MODULES} ${file_path} CACHE INTERNAL "GUDHI_MODULES")
  else()
    set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} ${file_path} CACHE INTERNAL "GUDHI_MISSING_MODULES")
  endif()
  # Required by user_version
  set(GUDHI_MODULES_FULL_LIST ${GUDHI_MODULES_FULL_LIST} ${file_path} PARENT_SCOPE)
  # Include module headers is independant - You may ask for no Alpha complex module but Python interface i.e.
  if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src/${file_path}/include/)
    include_directories(src/${file_path}/include/)
  endif()

endfunction(add_gudhi_module)

option(WITH_GUDHI_BENCHMARK "Activate/desactivate benchmark compilation" OFF)
option(WITH_GUDHI_EXAMPLE "Activate/desactivate examples compilation and installation" OFF)
option(WITH_GUDHI_PYTHON "Activate/desactivate python module compilation and installation" ON)
option(WITH_GUDHI_TEST "Activate/desactivate examples compilation and installation" ON)
option(WITH_GUDHI_UTILITIES "Activate/desactivate utilities compilation and installation" ON)

if (WITH_GUDHI_BENCHMARK)
  set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};benchmark")
endif()
if (WITH_GUDHI_EXAMPLE)
  set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};example")
endif()
if (WITH_GUDHI_TEST)
  set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};test")
endif()
if (WITH_GUDHI_UTILITIES)
  set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};utilities")
endif()

message("++ GUDHI_SUB_DIRECTORIES list is:\"${GUDHI_SUB_DIRECTORIES}\"")



