include_directories(${PROJECT_SOURCE_DIR}/src/)
set(SOURCE_FILES
  main.cpp
)

add_executable(read_stat ${SOURCE_FILES})
target_link_libraries(read_stat easehtslib)
