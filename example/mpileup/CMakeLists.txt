set(MPILEUP_FILES
        bam2bcf.c
        sample.c
        bedidx.c
        sam_opts.c
        sam_utils.c
        bam2bcf_indel.c
        kthread.c
        mcall.c
        vcmp.c
        ploidy.c
        regidx.c
)

add_executable(mpileup ${MPILEUP_FILES} mpileup.cpp)
add_executable(mpileup2 ${MPILEUP_FILES} mpileup2.cpp)
add_executable(mpileup_bcftools ${MPILEUP_FILES} mpileup_bcftools.cpp)
add_executable(mpileup_bcftools2 ${MPILEUP_FILES} mpileup_bcftools2.cpp)
add_executable(mpileup_bcftools_region ${MPILEUP_FILES} mpileup_bcftools_region.cpp)

include_directories(${PROJECT_SOURCE_DIR}/src/)
target_link_libraries(mpileup easehtslib)
target_link_libraries(mpileup2 easehtslib)
target_link_libraries(mpileup_bcftools easehtslib)
target_link_libraries(mpileup_bcftools2 easehtslib)
target_link_libraries(mpileup_bcftools_region easehtslib)



