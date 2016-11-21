include_directories(${PROJECT_SOURCE_DIR}/src/)
set(MUTECT_FILES
  mutect.cpp
  main.cpp
)

add_executable(mutect ${MUTECT_FILES})
target_link_libraries(mutect easehtslib)
