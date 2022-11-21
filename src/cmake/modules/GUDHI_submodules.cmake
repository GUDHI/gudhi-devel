# For those who dislike bundled dependencies, this indicates where to find a preinstalled Hera.
set(HERA_INTERNAL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ext/hera/include)
set(HERA_INCLUDE_DIR ${HERA_INTERNAL_INCLUDE_DIR} CACHE PATH "Directory where one can find hera/{wasserstein.h,bottleneck.h}")
