include_directories(${PROJECT_SOURCE_DIR}/src/)
add_executable(mutect mutect.cpp)
target_link_libraries(mutect easehtslib)
