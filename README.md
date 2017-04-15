# cpp-htslib
A wrapper of htslib which rewrite by cpp

# Building cpp-htslib

`sh build.sh` to build the project,and `sh build.sh test` to build and run unittest.
When clone the project, you should change the data dir of the test data

# NOTE
When you test the correctness of mutect, you should only set jdk1.7.0_79,
jdk1.7.0_80 is also not correct, Because the GATK framework use random to
produce the result. The mutect implentation also implentation the random
function as used in GATK, to make sure the answer is the same
