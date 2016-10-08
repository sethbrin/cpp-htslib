//
// Created by zp on 8/22/16.
//

#include "htslib/sam.h"

#include <string.h>
#include <gtest/gtest.h>

int fun(int a, int b) {
  return a + b;
}

TEST(fun_Test, TESTNORMAL) {
  EXPECT_EQ(fun(1, 2), 3);
}


int sam_read() {
  char filename[] = "/home/zp/work/rnaseq/cpp-htslib/test/uncompressed.sam";
  samFile *in = sam_open(filename, "r");
  if (in == NULL) {
    fprintf(stderr, "Error opening \"%s\"\n", filename);
    return EXIT_FAILURE;
  }

  bam_hdr_t *h = sam_hdr_read(in);
  if (h == NULL) {
    fprintf(stderr, "Couldn't read header for \"%s\"\n", filename);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

TEST(sam_read, header) {
  EXPECT_EQ(sam_read(), 0);

}
