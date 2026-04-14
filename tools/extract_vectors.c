// Tool to extract intermediate test vectors from the C implementation
// for cross-checking the Rust port. Links against the lvl1 build.
// Build: see tools/build_extract.sh

#include <stdio.h>
#include <string.h>
#include <fp.h>
#include <fp2.h>
#include <rng.h>

static void hex(const char *name, const void *p, size_t n) {
    const unsigned char *b = p;
    printf("%s = ", name);
    for (size_t i = 0; i < n; i++) printf("%02x", b[i]);
    printf("\n");
}

int main(void) {
    // DRBG vectors
    unsigned char seed[48] = {0};
    for (int i = 0; i < 48; i++) seed[i] = (unsigned char)i;
    randombytes_init(seed, NULL, 256);
    unsigned char buf[64];
    randombytes(buf, 64);
    hex("drbg_seq_seed_out64", buf, 64);
    randombytes(buf, 32);
    hex("drbg_seq_seed_out32", buf, 32);

    // Fp encode/decode round-trip
    fp_t a, b, c;
    unsigned char enc[FP_ENCODED_BYTES];
    fp_set_one(&a);
    fp_encode(enc, &a);
    hex("fp_one_encoded", enc, FP_ENCODED_BYTES);

    fp_set_small(&a, 5);
    fp_set_small(&b, 7);
    fp_mul(&c, &a, &b);
    fp_encode(enc, &c);
    hex("fp_5x7_encoded", enc, FP_ENCODED_BYTES);

    fp_inv(&c);
    fp_encode(enc, &c);
    hex("fp_5x7_inv_encoded", enc, FP_ENCODED_BYTES);

    // Fp2: (3+4i)^2 = -7 + 24i
    fp2_t x, y;
    unsigned char enc2[FP2_ENCODED_BYTES];
    fp_set_small(&x.re, 3);
    fp_set_small(&x.im, 4);
    fp2_sqr(&y, &x);
    fp2_encode(enc2, &y);
    hex("fp2_3_4i_sq_encoded", enc2, FP2_ENCODED_BYTES);

    return 0;
}
