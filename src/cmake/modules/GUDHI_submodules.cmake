# For those who dislike bundled dependencies, this indicates where to find a preinstalled Hera.
set(HERA_INTERNAL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ext/hera/include)
set(HERA_INCLUDE_DIR ${HERA_INTERNAL_INCLUDE_DIR} CACHE PATH "Directory where one can find hera/{wasserstein.h,bottleneck.h}")
# since everything is cleanly under include/hera/, there is no harm always including it
include_directories(${HERA_INCLUDE_DIR})

if (NOT EXISTS ${HERA_INCLUDE_DIR}/hera/wasserstein.h OR NOT EXISTS ${HERA_INCLUDE_DIR}/hera/bottleneck.h)
  message(WARNING "${HERA_INCLUDE_DIR}/hera/{wasserstein.h,bottleneck.h} are not found.\n\
  GUDHI requires this submodules, please consider to launch `git submodule update --init`.\n\
  If hera was installed in a specific directory, you can also consider to specify it to the cmake command with `cmake -DHERA_INCLUDE_DIR=... ...`")
endif()