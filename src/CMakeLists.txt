set(SOURCE_FILES
    utils.cpp
    sam_bam_reader.cpp
    sam_bam_record.cpp
    genome_loc.cpp
    reference_sequence.cpp
)

add_library(easehtslib ${SOURCE_FILES})
target_link_libraries(easehtslib ${HTS_LIB}/libhts.a)

add_executable(main main.cpp)
link_directories(${LIBRARY_OUTPUT_PATH})
target_link_libraries(main easehtslib)

if (test)
    link_libraries(gtest gtest_main)

    create_test(main_unittest)
    create_test(sam_bam_reader_unittest)
    create_test(genome_loc_unittest)
    create_test(reference_sequence_unittest)
    create_test(pileup_unittest)
endif()
