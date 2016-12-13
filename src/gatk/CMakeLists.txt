if (test)
    link_libraries(gtest gtest_main)

    create_test(alignment_state_machine_unittest)
endif()
