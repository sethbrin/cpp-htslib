#ifndef UNIT_TEST_H_
#define UNIT_TEST_H_

#include <string.h>

char test_dir[] = "/home/zp/work/rnaseq/cpp-htslib/test/";

// 注意此处定义var，TEST_FILE只能单独作为一个语句,
// 嵌套进if等块中可能会出问题
#define TEST_FILE(filename, var) \
  char var[1024]; \
  do { \
    strcpy(var, test_dir); \
    strcat(var, filename); \
  } while(0)

#endif
