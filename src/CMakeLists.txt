set(SOURCE_FILES
    utils.cpp
    sam_seqence_record.cpp
    sam_bam_normal_reader.cpp
    reader.cpp
    sam_bam_record.cpp
    bam_index_reader.cpp
    bam_index.cpp
    genome_loc.cpp
    reference_sequence.cpp
)

add_library(easehtslib ${SOURCE_FILES})
target_link_libraries(easehtslib ${HTS_LIB}/libhts.a)

add_executable(main main.cpp)
link_directories(${LIBRARY_OUTPUT_PATH})
target_link_libraries(main easehtslib)

set(MPILEUP_FILES
        mpileup2.cpp
        bam2bcf.c
        sample.c
        bedidx.c
        sam_opts.c
        sam_utils.c
        pileup_buffer.cpp
        bam2bcf_indel.c
)
add_executable(mpileup ${MPILEUP_FILES})
target_link_libraries(mpileup easehtslib)


if (test)
    link_libraries(gtest gtest_main)

    create_test(main_unittest)
    create_test(reader_unittest)
    create_test(sam_bam_read_unittest)
    create_test(sam_bam_record_unittest)
    create_test(genome_loc_unittest)
    create_test(reference_sequence_unittest)
endif()
