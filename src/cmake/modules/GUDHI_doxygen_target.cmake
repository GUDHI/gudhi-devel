# add a target to generate API documentation with Doxygen
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
  set(GUDHI_MODULES ${GUDHI_MODULES} "cpp-documentation" CACHE INTERNAL "GUDHI_MODULES")

  # starting from cmake 3.9 the usage of DOXYGEN_EXECUTABLE is deprecated
  if(TARGET Doxygen::doxygen)
    get_property(DOXYGEN_EXECUTABLE TARGET Doxygen::doxygen PROPERTY IMPORTED_LOCATION)
  endif()

  add_custom_target(doxygen ${DOXYGEN_EXECUTABLE} ${GUDHI_USER_VERSION_DIR}/Doxyfile
                    WORKING_DIRECTORY ${GUDHI_USER_VERSION_DIR}
                    COMMENT "Generating API documentation with Doxygen in ${GUDHI_USER_VERSION_DIR}/doc/html/" VERBATIM)

  if(TARGET user_version)
    # In dev version, doxygen target depends on user_version target. Not existing in user version
    add_dependencies(doxygen user_version)
  endif()
else()
  set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} "cpp-documentation" CACHE INTERNAL "GUDHI_MISSING_MODULES")
endif()
