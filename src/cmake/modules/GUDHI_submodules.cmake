# For those who dislike bundled dependencies, this indicates where to find a preinstalled Hera.
set(HERA_WASSERSTEIN_INTERNAL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ext/hera/wasserstein/include)
set(HERA_WASSERSTEIN_INCLUDE_DIR ${HERA_WASSERSTEIN_INTERNAL_INCLUDE_DIR} CACHE PATH "Directory where one can find Hera's wasserstein.h")
set(HERA_BOTTLENECK_INTERNAL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ext/hera/bottleneck/include)
set(HERA_BOTTLENECK_INCLUDE_DIR ${HERA_BOTTLENECK_INTERNAL_INCLUDE_DIR} CACHE PATH "Directory where one can find Hera's bottleneck.h")