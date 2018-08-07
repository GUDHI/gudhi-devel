# Definition of the custom target user_version
add_custom_target(user_version)

if(DEFINED USER_VERSION_DIR)
  # set the GUDHI_USER_VERSION_DIR with USER_VERSION_DIR defined by the user
  set(GUDHI_USER_VERSION_DIR ${CMAKE_CURRENT_BINARY_DIR}/${USER_VERSION_DIR})
else()
  # set the GUDHI_USER_VERSION_DIR with timestamp and Gudhi version number
  string(TIMESTAMP DATE_AND_TIME "%Y-%m-%d-%H-%M-%S")
  set(GUDHI_USER_VERSION_DIR ${CMAKE_CURRENT_BINARY_DIR}/${DATE_AND_TIME}_GUDHI_${GUDHI_VERSION})
endif()

add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   make_directory ${GUDHI_USER_VERSION_DIR}
                   COMMENT "user_version creation in ${GUDHI_USER_VERSION_DIR}")

add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy ${CMAKE_SOURCE_DIR}/Conventions.txt ${GUDHI_USER_VERSION_DIR}/Conventions.txt)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy ${CMAKE_SOURCE_DIR}/README ${GUDHI_USER_VERSION_DIR}/README)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy ${CMAKE_SOURCE_DIR}/COPYING ${GUDHI_USER_VERSION_DIR}/COPYING)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy ${CMAKE_SOURCE_DIR}/src/CMakeLists.txt ${GUDHI_USER_VERSION_DIR}/CMakeLists.txt)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy ${CMAKE_SOURCE_DIR}/src/Doxyfile ${GUDHI_USER_VERSION_DIR}/Doxyfile)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy ${CMAKE_SOURCE_DIR}/src/GUDHIConfigVersion.cmake.in ${GUDHI_USER_VERSION_DIR}/GUDHIConfigVersion.cmake.in)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy ${CMAKE_SOURCE_DIR}/src/GUDHIConfig.cmake.in ${GUDHI_USER_VERSION_DIR}/GUDHIConfig.cmake.in)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy ${CMAKE_SOURCE_DIR}/CMakeGUDHIVersion.txt ${GUDHI_USER_VERSION_DIR}/CMakeGUDHIVersion.txt)

add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy_directory ${CMAKE_SOURCE_DIR}/biblio ${GUDHI_USER_VERSION_DIR}/biblio)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy_directory ${CMAKE_SOURCE_DIR}/src/cython ${GUDHI_USER_VERSION_DIR}/cython)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy_directory ${CMAKE_SOURCE_DIR}/data ${GUDHI_USER_VERSION_DIR}/data)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy_directory ${CMAKE_SOURCE_DIR}/src/cmake ${GUDHI_USER_VERSION_DIR}/cmake)
add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                   copy_directory ${CMAKE_SOURCE_DIR}/src/GudhUI ${GUDHI_USER_VERSION_DIR}/GudhUI)

set(GUDHI_DIRECTORIES "doc;example;concept;utilities")
if (CGAL_VERSION VERSION_LESS 4.11.0)
  set(GUDHI_INCLUDE_DIRECTORIES "include/gudhi;include/gudhi_patches")
else ()
    set(GUDHI_INCLUDE_DIRECTORIES "include/gudhi")
endif ()

foreach(GUDHI_MODULE ${GUDHI_MODULES_FULL_LIST})
  foreach(GUDHI_DIRECTORY ${GUDHI_DIRECTORIES})
    # Find files
    file(GLOB GUDHI_FILES ${CMAKE_SOURCE_DIR}/src/${GUDHI_MODULE}/${GUDHI_DIRECTORY}/*)

    foreach(GUDHI_FILE ${GUDHI_FILES})
      get_filename_component(GUDHI_FILE_NAME ${GUDHI_FILE} NAME)
      # GUDHI_FILE can be a file or a directory
      if(IS_DIRECTORY ${GUDHI_FILE})
        add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                           copy_directory ${GUDHI_FILE} ${GUDHI_USER_VERSION_DIR}/${GUDHI_DIRECTORY}/${GUDHI_MODULE}/${GUDHI_FILE_NAME})
      else()
        add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                           copy ${GUDHI_FILE} ${GUDHI_USER_VERSION_DIR}/${GUDHI_DIRECTORY}/${GUDHI_MODULE}/${GUDHI_FILE_NAME})
      endif()
    endforeach()
  endforeach(GUDHI_DIRECTORY ${GUDHI_DIRECTORIES})

  foreach(GUDHI_INCLUDE_DIRECTORY ${GUDHI_INCLUDE_DIRECTORIES})
    # include files
    file(GLOB GUDHI_INCLUDE_FILES ${CMAKE_SOURCE_DIR}/src/${GUDHI_MODULE}/${GUDHI_INCLUDE_DIRECTORY}/*)

    foreach(GUDHI_INCLUDE_FILE ${GUDHI_INCLUDE_FILES})
      get_filename_component(GUDHI_INCLUDE_FILE_NAME ${GUDHI_INCLUDE_FILE} NAME)
      # GUDHI_INCLUDE_FILE can be a file or a directory
      if(IS_DIRECTORY ${GUDHI_INCLUDE_FILE})
        add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                           copy_directory ${GUDHI_INCLUDE_FILE} ${GUDHI_USER_VERSION_DIR}/${GUDHI_INCLUDE_DIRECTORY}/${GUDHI_INCLUDE_FILE_NAME})
      else()
        add_custom_command(TARGET user_version PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
                           copy ${GUDHI_INCLUDE_FILE} ${GUDHI_USER_VERSION_DIR}/${GUDHI_INCLUDE_DIRECTORY}/${GUDHI_INCLUDE_FILE_NAME})
      endif()
    endforeach()
  endforeach(GUDHI_INCLUDE_DIRECTORY ${GUDHI_INCLUDE_DIRECTORIES})

endforeach(GUDHI_MODULE ${GUDHI_MODULES_FULL_LIST})