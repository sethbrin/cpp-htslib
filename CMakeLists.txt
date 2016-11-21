cmake_minimum_required(VERSION 2.6)
project(easehtslib)

#set(CMAKE_CXX_COMPILER clang++)

set(CXX_FLAGS
        -g
        # -DVALGRIND
        # -DMUDUO_STD_STRING
        -DCHECK_PTHREAD_RETURN_VALUE
        -D_FILE_OFFSET_BITS=64
        -fpermissive
        #-Wall
        #-Wextra
        #-Werror
        #-Wno-unused-parameter
        #-Woverloaded-virtual
        #-Wpointer-arith
        #-Wwrite-strings
        #-march=native
        # -MMD
        #-std=c++0x
        #-fno-elide-constructors
        -std=c++11
        -rdynamic
        )
# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(test "Build all tests." ON) # Makes boolean 'test' available.
option(EASEHTS_BUILD_THIRDPARTY_GTEST
  "Use gtest installation in `thirdparty/gtest` by default if available" ON)

set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -finline-limit=1000 -DNDEBUG")
string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)


find_package(ZLIB)

set(HTS_LIB "${PROJECT_SOURCE_DIR}/thirdparty/htslib/")

include_directories(${HTS_LIB})

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
