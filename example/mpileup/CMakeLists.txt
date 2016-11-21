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
include_directories(${PROJECT_SOURCE_DIR}/src/)
target_link_libraries(mpileup easehtslib)



