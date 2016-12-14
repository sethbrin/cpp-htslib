SET(BENCHMARK_SEARCH_PATH
  "${BENCHMARK_SOURCE_DIR}"
  "${CMAKE_CURRENT_LIST_DIR}/../thirdparty/benchmark")

FIND_PATH(GTEST_SOURCE_DIR
    NAMES CMakeLists.txt src/benchmark.cc
    PATHS ${BENCHMARK_SEARCH_PATH})


# Debian installs gtest include directory in /usr/include, thus need to look
# for include directory separately from source directory.
FIND_PATH(BENCHMARK_INCLUDE_DIR
    NAMES include/benchmark.h
    PATH_SUFFIXES include
    HINTS ${BENCHMARK_SOURCE_DIR}
    PATHS ${BENCHMARK_SEARCH_PATH})

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BenchmarkSrc DEFAULT_MSG
  BENCHMARK_SOURCE_DIR
  BENCHMARK_INCLUDE_DIR)
