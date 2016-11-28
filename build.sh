#!/bin/sh

#set -x

SOURCE_DIR=`pwd`
BUILD_DIR=${BUILD_DIR:-./build}
BUILD_TYPE=${BUILD_TYPE:-./release}
#BUILD_TYPE=${BUILD_TYPE:-./debug}
INSTALL_DIR=${INSTALL_DIR:-${BUILD_TYPE}-install}

# build htslib
cd htslib \
  && make \
  && cd ..

mkdir -p $BUILD_DIR/$BUILD_TYPE \
  && cd $BUILD_DIR/$BUILD_TYPE \
  && cmake \
           -Dtest=ON \
           -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
           -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
           $SOURCE_DIR \
  && make $*

# 生成单元测试头文件，指定测试文件路径
cat <<EOF > ${SOURCE_DIR}"/src/easehts/unittest.h"
#ifndef UNIT_TEST_H_
#define UNIT_TEST_H_

#include <string.h>

char test_dir[] = "${SOURCE_DIR}/test/";

// 注意此处定义var，TEST_FILE只能单独作为一个语句,
// 嵌套进if等块中可能会出问题
#define TEST_FILE(filename, var) \\
  char var[1024]; \\
  do { \\
    strcpy(var, test_dir); \\
    strcat(var, filename); \\
  } while(0)

#endif
EOF

if [ $# -gt 0 ] && [ $1 == "test" ];then
  for unittest in bin/*_unittest; do
    echo "=======run ${unittest}==========="
    ./${unittest}
    echo && echo
  done
fi
