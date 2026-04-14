#!/bin/bash
set -e
C=/root/src/personal-hacking/the-sqisign
B=$C/build
gcc -O2 -o /tmp/dump_keygen \
  /root/src/personal-hacking/sqisign-rs/tools/dump_keygen.c \
  -I$C/include \
  -I$C/src/gf/ref/include \
  -I$C/src/ec/ref/include \
  -I$C/src/hd/ref/include \
  -I$C/src/mp/ref/generic \
  -I$C/src/precomp/ref/lvl1/include \
  -I$C/src/common/generic/include \
  -I$C/src/common/ref/include \
  -I$C/src/verification/ref/include \
  -I$C/src/signature/ref/include \
  -I$C/src/quaternion/ref/generic/include \
  -I$C/src/id2iso/ref/include \
  -I$C/src/mini-gmp \
  -DSQISIGN_VARIANT=lvl1 -DRADIX_64 -DNDEBUG -DSQISIGN_BUILD_TYPE_REF -DSQISIGN_GF_IMPL_REF -DENABLE_SIGN -DMINI_GMP \
  -Wl,--start-group \
  $B/src/libsqisign_lvl1_test.a \
  $B/src/signature/ref/lvl1/libsqisign_signature_lvl1.a \
  $B/src/id2iso/ref/lvl1/libsqisign_id2iso_lvl1.a \
  $B/src/verification/ref/lvl1/libsqisign_verification_lvl1.a \
  $B/src/hd/ref/lvl1/libsqisign_hd_lvl1.a \
  $B/src/ec/ref/lvl1/libsqisign_ec_lvl1.a \
  $B/src/gf/ref/lvl1/libsqisign_gf_lvl1.a \
  $B/src/mp/ref/generic/libsqisign_mp_generic.a \
  $B/src/precomp/ref/lvl1/libsqisign_precomp_lvl1.a \
  $B/src/quaternion/ref/generic/libsqisign_quaternion_generic.a \
  $B/src/common/generic/libsqisign_common_test.a \
  $B/src/common/generic/libsqisign_common_sys.a \
  $B/libGMP.a \
  -Wl,--end-group -lm
/tmp/dump_keygen
