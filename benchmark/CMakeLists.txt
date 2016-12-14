include_directories(${PROJECT_SOURCE_DIR}/src/)

find_package(BenchmarkSrc)
include_directories(SYSTEM ${BENCHMARK_INCLUDE_DIR})

link_libraries(easehtslib)
link_libraries(benchmark)

add_executable(unordered_map_bench unordered_map_bench.cpp)
add_executable(inline_bench inline_bench.cpp)

