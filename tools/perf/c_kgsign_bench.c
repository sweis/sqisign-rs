// Controlled C benchmark for keygen/sign on KAT-0 seed.
// Modes (argv[1]): "kg"     -> 20x { reseed; keypair }
//                  "kgsign" -> 20x { reseed; keypair; sign(KAT0_MSG) }
// Each iteration is identical (DRBG re-seeded), so perf instruction counts
// are directly comparable to the Rust harness doing the same.

#include <stdio.h>
#include <string.h>
#include <api.h>
#include <rng.h>

static unsigned char SEED[48] = {
  0x06,0x15,0x50,0x23,0x4D,0x15,0x8C,0x5E,0xC9,0x55,0x95,0xFE,0x04,0xEF,0x7A,0x25,
  0x76,0x7F,0x2E,0x24,0xCC,0x2B,0xC4,0x79,0xD0,0x9D,0x86,0xDC,0x9A,0xBC,0xFD,0xE7,
  0x05,0x6A,0x8C,0x26,0x6F,0x9E,0xF9,0x7E,0xD0,0x85,0x41,0xDB,0xD2,0xE1,0xFF,0xA1,
};
static const unsigned char MSG[33] = {
  0xD8,0x1C,0x4D,0x8D,0x73,0x4F,0xCB,0xFB,0xEA,0xDE,0x3D,0x3F,0x8A,0x03,0x9F,0xAA,
  0x2A,0x2C,0x99,0x57,0xE8,0x35,0xAD,0x55,0xB2,0x2E,0x75,0xBF,0x57,0xBB,0x55,0x6A,
  0xC8,
};

#define ITERS 20

int main(int argc, char **argv) {
    int do_sign = (argc > 1 && strcmp(argv[1], "kgsign") == 0);
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char sm[CRYPTO_BYTES + sizeof(MSG)];
    unsigned long long smlen;

    for (int i = 0; i < ITERS; i++) {
        randombytes_init(SEED, NULL, 256);
        if (crypto_sign_keypair(pk, sk) != 0) { fprintf(stderr, "kp fail\n"); return 1; }
        if (do_sign) {
            if (crypto_sign(sm, &smlen, MSG, sizeof(MSG), sk) != 0) { fprintf(stderr, "sign fail\n"); return 1; }
        }
    }
    // Prevent optimizer from eliding work.
    volatile unsigned char sink = pk[0] ^ sk[0] ^ (do_sign ? sm[0] : 0);
    (void)sink;
    return 0;
}
