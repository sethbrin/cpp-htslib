#!/bin/sh

#set -x

SOURCE_DIR=`pwd`
BUILD_DIR=${BUILD_DIR:-./build}
#BUILD_TYPE=${BUILD_TYPE:-./release}
BUILD_TYPE=${BUILD_TYPE:-./debug}
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

if [ $# -gt 0 ] && [ $1 == "test" ];then
  for unittest in bin/*_unittest; do
    echo "=======run ${unittest}==========="
    ./${unittest}
    echo && echo
  done
fi
