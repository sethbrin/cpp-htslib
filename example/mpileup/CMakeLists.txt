set(MPILEUP_FILES
        bam2bcf.c
        sample.c
        bedidx.c
        sam_opts.c
        sam_utils.c
        bam2bcf_indel.c
        kthread.c
)

add_executable(mpileup ${MPILEUP_FILES} mpileup.cpp)
add_executable(mpileup2 ${MPILEUP_FILES} mpileup2.cpp)

include_directories(${PROJECT_SOURCE_DIR}/src/)
target_link_libraries(mpileup easehtslib)
target_link_libraries(mpileup2 easehtslib)



