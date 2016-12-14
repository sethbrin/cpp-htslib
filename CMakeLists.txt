cmake_minimum_required(VERSION 2.6)
project(easehtslib)

#set(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries_2016/linux/bin/intel64/icpc)

set(CXX_FLAGS
        -g
        -O3
        #-pg
        -fpermissive
        #-Wall
        #-Wextra
        #-Werror
        -march=native
        -finline-functions
        -Wno-unused-parameter
        -Woverloaded-virtual
        -Wpointer-arith
        -Wwrite-strings
        -fno-elide-constructors
        -std=c++11
        -rdynamic
        )
# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(test "Build all tests." ON) # Makes boolean 'test' available.
option(EASEHTS_BUILD_THIRDPARTY_GTEST
  "Use gtest installation in `thirdparty/gtest` by default if available" ON)
option(USE_TCMALLOC "use Tcmalloc." ON)
option(BUILD_BENCHMARK "build benchmark example." ON)

set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -finline-limit=1000 -DNDEBUG")
string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

if (USE_TCMALLOC)
  link_libraries(tcmalloc)
endif()

#find_package(ZLIB)

set(HTS_LIB "${PROJECT_SOURCE_DIR}/thirdparty/htslib/")
include_directories(${HTS_LIB})

find_package(Tcmalloc)


link_libraries(
        pthread
        z
        unwind
)

if (test)
    macro (create_test target)
        add_executable(${target} ${target}.cpp)
        target_link_libraries(${target} easehtslib)
        add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}>)
    endmacro (create_test)

    enable_testing()

    find_package(GTestSrc)
    add_subdirectory(${GTEST_SOURCE_DIR} ${CMAKE_BINARY_DIR}/googletest)
    include_directories(SYSTEM ${GTEST_INCLUDE_DIR})

endif()

add_subdirectory(src)
add_subdirectory(example)

if (BUILD_BENCHMARK)
  add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/benchmark)
  add_subdirectory(benchmark)
endif()
