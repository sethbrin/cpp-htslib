include_directories(${PROJECT_SOURCE_DIR}/src/)
set(SOURCE_FILES
  main.cpp
)

add_executable(coverage ${SOURCE_FILES})
target_link_libraries(coverage easehtslib)
