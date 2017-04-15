include_directories(${PROJECT_SOURCE_DIR}/src/)
set(SOURCE_FILES
  main.cpp
)

add_executable(merge_interval ${SOURCE_FILES})
target_link_libraries(merge_interval easehtslib)
