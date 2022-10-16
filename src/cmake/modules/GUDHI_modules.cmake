# A function to add a new module in GUDHI

set(GUDHI_MODULES_FULL_LIST "")
function(add_gudhi_module file_path)
  option("WITH_MODULE_GUDHI_${file_path}" "Activate/deactivate ${file_path} compilation and installation" ON)
  if (WITH_MODULE_GUDHI_${file_path})
    set(GUDHI_MODULES ${GUDHI_MODULES} ${file_path} CACHE INTERNAL "GUDHI_MODULES")
  else()
    set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} ${file_path} CACHE INTERNAL "GUDHI_MISSING_MODULES")
  endif()
  # Required by user_version
  set(GUDHI_MODULES_FULL_LIST ${GUDHI_MODULES_FULL_LIST} ${file_path} PARENT_SCOPE)
  # Include module headers is independent - You may ask for no Alpha complex module but Python interface i.e.
  if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src/${file_path}/include/)
    include_directories(src/${file_path}/include/)
  endif()

endfunction(add_gudhi_module)

if (WITH_GUDHI_BENCHMARK)
  set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};benchmark")
endif()
if (WITH_GUDHI_EXAMPLE)
  set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};example")
endif()
if (WITH_GUDHI_TEST)
  # All tests are using boost tests
  if(TARGET Boost::unit_test_framework)
    set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};test")
  else()
    message("++ WITH_GUDHI_TEST but no TARGET Boost::unit_test_framework")
  endif()
endif()
if (WITH_GUDHI_UTILITIES)
  set(GUDHI_SUB_DIRECTORIES "${GUDHI_SUB_DIRECTORIES};utilities")
endif()

message("++ GUDHI_SUB_DIRECTORIES list is:\"${GUDHI_SUB_DIRECTORIES}\"")



