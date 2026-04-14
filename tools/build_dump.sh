#!/bin/bash
set -e
C=/root/src/personal-hacking/the-sqisign
B=$C/build
gcc -O2 -o /tmp/dump_verify \
  /root/src/personal-hacking/sqisign-rs/tools/dump_verify_intermediates.c \
  -I$C/include \
  -I$C/src/gf/ref/include \
  -I$C/src/ec/ref/include \
  -I$C/src/hd/ref/include \
  -I$C/src/mp/ref/generic \
  -I$C/src/precomp/ref/lvl1/include \
  -I$C/src/common/generic/include \
  -I$C/src/common/ref/include \
  -I$C/src/verification/ref/include \
  -DSQISIGN_VARIANT=lvl1 -DRADIX_64 -DNDEBUG -DSQISIGN_BUILD_TYPE_REF -DSQISIGN_GF_IMPL_REF \
  -Wl,--start-group \
  $B/src/verification/ref/lvl1/libsqisign_verification_lvl1.a \
  $B/src/hd/ref/lvl1/libsqisign_hd_lvl1.a \
  $B/src/ec/ref/lvl1/libsqisign_ec_lvl1.a \
  $B/src/gf/ref/lvl1/libsqisign_gf_lvl1.a \
  $B/src/mp/ref/generic/libsqisign_mp_generic.a \
  $B/src/precomp/ref/lvl1/libsqisign_precomp_lvl1.a \
  $B/src/common/generic/libsqisign_common_sys.a \
  -Wl,--end-group
/tmp/dump_verify
