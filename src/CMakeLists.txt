set(SOURCE_FILES
    args_parser.cpp
    genome_loc.cpp
    utils.cpp
    sam_bam_reader.cpp
    sam_bam_record.cpp
    reference_sequence.cpp
    pileup.cpp
    binomial_distribution.cpp
    base_utils.cpp
    diploid_SNP_genotype_likelihoods.cpp
    diploid_genotype.cpp
    gatk/alignment_state_machine.cpp
    gatk/pileup.cpp
)

add_library(easehtslib ${SOURCE_FILES})
target_link_libraries(easehtslib ${HTS_LIB}/libhts.a)
include_directories(${PROJECT_SOURCE_DIR}/src/)

add_subdirectory(gatk)

if (test)
    link_libraries(gtest gtest_main)

    create_test(sam_bam_reader_unittest)
    create_test(genome_loc_unittest)
    create_test(reference_sequence_unittest)
    create_test(pileup_unittest)
    create_test(sam_bam_record_unittest)
    create_test(utils_unittest)
endif()
