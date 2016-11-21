set(SOURCE_FILES
    args_parser.cpp
    genome_loc.cpp
    utils.cpp
    sam_bam_reader.cpp
    sam_bam_record.cpp
    reference_sequence.cpp
)

add_library(easehtslib ${SOURCE_FILES})
target_link_libraries(easehtslib ${HTS_LIB}/libhts.a)

if (test)
    link_libraries(gtest gtest_main)

    create_test(sam_bam_reader_unittest)
    create_test(genome_loc_unittest)
    create_test(reference_sequence_unittest)
    create_test(pileup_unittest)
    create_test(sam_bam_record_unittest)
endif()
