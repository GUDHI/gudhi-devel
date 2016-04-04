# Install script for directory: /home/frg/Bureau/mWorkingDirectory

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Simplex_tree/test/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Simplex_tree/example/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/test/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Persistent_cohomology/example/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Skeleton_blocker/test/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Skeleton_blocker/example/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Contraction/example/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Hasse_complex/example/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Alpha_shapes/example/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Alpha_shapes/test/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Bipartite_graphs_matching/example/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/src/Bipartite_graphs_matching/test/cmake_install.cmake")
  include("/home/frg/Bureau/mWorkingDirectory/build/data/points/generator/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

file(WRITE "/home/frg/Bureau/mWorkingDirectory/build/${CMAKE_INSTALL_MANIFEST}" "")
foreach(file ${CMAKE_INSTALL_MANIFEST_FILES})
  file(APPEND "/home/frg/Bureau/mWorkingDirectory/build/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
endforeach()
