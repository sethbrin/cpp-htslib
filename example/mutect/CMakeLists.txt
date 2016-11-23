include_directories(${PROJECT_SOURCE_DIR}/src/)
set(MUTECT_FILES
  mutect.cpp
  locus_read_pile.cpp
  main.cpp
)

add_executable(mutect ${MUTECT_FILES})
target_link_libraries(mutect easehtslib)

if (test)
    link_libraries(gtest gtest_main)

    create_test(power_calculator_unittest)
endif()
