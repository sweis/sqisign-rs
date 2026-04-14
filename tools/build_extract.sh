#!/bin/bash
set -e
C=/root/src/personal-hacking/the-sqisign
B=$C/build
gcc -O2 -o /tmp/extract_vectors \
  /root/src/personal-hacking/sqisign-rs/tools/extract_vectors.c \
  -I$C/include \
  -I$C/src/gf/ref/include \
  -I$C/src/precomp/ref/lvl1/include \
  -I$C/src/common/generic/include \
  -I$C/src/common/ref/include \
  -DSQISIGN_VARIANT=1 -DRADIX_64 -DENABLE_SIGN \
  $B/src/libsqisign_lvl1.a \
  $B/src/libsqisign_gf_lvl1.a \
  $B/src/libsqisign_common_sys.a \
  $B/src/libsqisign_common_test.a \
  2>&1 || find $B -name "*.a" | head -20
/tmp/extract_vectors
