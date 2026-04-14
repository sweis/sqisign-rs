// Dump intermediate DRBG byte-counts and ideal norms during keygen for KAT 0.
#include <stdio.h>
#include <string.h>
#include <signature.h>
#include <quaternion_data.h>
#include <quaternion_constants.h>
#include <torsion_constants.h>
#include <id2iso.h>
#include <rng.h>
#include <stdlib.h>
static void pz(const char *name, const mpz_t z) {
    char *s = mpz_get_str(NULL, 10, z);
    printf("%s = %s\n", name, s);
    free(s);
}
#define gmp_printf(fmt, z) pz(fmt, z)

static unsigned long g_bytes = 0;
extern int randombytes(unsigned char *buf, unsigned long long n);
// Wrap randombytes to count calls — done by interposition below.

static void hex(const char *name, const void *p, size_t n) {
    const unsigned char *b = p;
    printf("%s = ", name);
    for (size_t i = 0; i < n; i++) printf("%02x", b[i]);
    printf("\n");
}

int main(void) {
    unsigned char seed[48] = {
        0x06,0x15,0x50,0x23,0x4D,0x15,0x8C,0x5E,0xC9,0x55,0x95,0xFE,0x04,0xEF,0x7A,0x25,
        0x76,0x7F,0x2E,0x24,0xCC,0x2B,0xC4,0x79,0xD0,0x9D,0x86,0xDC,0x9A,0xBC,0xFD,0xE7,
        0x05,0x6A,0x8C,0x26,0x6F,0x9E,0xF9,0x7E,0xD0,0x85,0x41,0xDB,0xD2,0xE1,0xFF,0xA1
    };
    randombytes_init(seed, NULL, 256);

    secret_key_t sk;
    public_key_t pk = {0};
    secret_key_init(&sk);

    printf("=== keygen ===\n");
    // Replicate the body of quat_sampling_random_ideal_O0_given_norm step by step
    {
        ibz_t n_temp, norm_d, disc;
        quat_alg_elem_t gen, gen_rerand;
        ibz_init(&n_temp); ibz_init(&norm_d); ibz_init(&disc);
        quat_alg_elem_init(&gen); quat_alg_elem_init(&gen_rerand);
        int iter = 0;
        int found = 0;
        while (!found) {
            iter++;
            ibz_set(&gen.coord[0], 0);
            ibz_sub(&n_temp, &SEC_DEGREE, &ibz_const_one);
            for (int i = 1; i < 4; i++)
                ibz_rand_interval(&gen.coord[i], &ibz_const_zero, &n_temp);
            quat_alg_norm(&n_temp, &norm_d, &gen, &QUATALG_PINFTY);
            ibz_neg(&disc, &n_temp);
            ibz_mod(&disc, &disc, &SEC_DEGREE);
            found = ibz_sqrt_mod_p(&gen.coord[0], &disc, &SEC_DEGREE);
            found = found && !quat_alg_elem_is_zero(&gen);
        }
        printf("phase1_iters = %d\n", iter);
        pz("gen.coord[0]", gen.coord[0]);
        pz("gen.coord[1]", gen.coord[1]);
        pz("gen.coord[2]", gen.coord[2]);
        pz("gen.coord[3]", gen.coord[3]);
        unsigned char r[8]; randombytes(r, 8);
        hex("drbg_after_phase1", r, 8);
        ibz_finalize(&n_temp); ibz_finalize(&norm_d); ibz_finalize(&disc);
        quat_alg_elem_finalize(&gen); quat_alg_elem_finalize(&gen_rerand);
    }

    randombytes_init(seed, NULL, 256);
    int found = quat_sampling_random_ideal_O0_given_norm(
        &sk.secret_ideal, &SEC_DEGREE, 1, &QUAT_represent_integer_params, NULL);
    printf("after_sampling_found = %d\n", found);
    pz("after_sampling_norm", sk.secret_ideal.norm);
    pz("after_sampling_basis_00", sk.secret_ideal.lattice.basis[0][0]);
    pz("after_sampling_basis_33", sk.secret_ideal.lattice.basis[3][3]);
    pz("after_sampling_denom", sk.secret_ideal.lattice.denom);

    unsigned char r[8];
    randombytes(r, 8);
    hex("drbg_after_sampling_next8", r, 8);

    // Re-seed and do full keygen for comparison
    randombytes_init(seed, NULL, 256);
    secret_key_finalize(&sk);
    secret_key_init(&sk);

    found = quat_sampling_random_ideal_O0_given_norm(
        &sk.secret_ideal, &SEC_DEGREE, 1, &QUAT_represent_integer_params, NULL);
    // Replicate prime_norm_reduced_equivalent body
    {
        ibz_mat_4x4_t gram, red, gin;
        ibz_mat_4x4_init(&gram); ibz_mat_4x4_init(&red); ibz_mat_4x4_init(&gin);
        quat_lideal_class_gram(&gin, &sk.secret_ideal, &QUATALG_PINFTY);
        for (int i=0;i<4;i++) for (int j=0;j<4;j++) {
            char nm[32]; snprintf(nm,32,"gin[%d][%d]",i,j);
            pz(nm, gin[i][j]);
        }
        quat_lideal_reduce_basis(&red, &gram, &sk.secret_ideal, &QUATALG_PINFTY);
        for (int i=0;i<4;i++) for (int j=0;j<4;j++) {
            char nm[32]; snprintf(nm,32,"red[%d][%d]",i,j);
            pz(nm, red[i][j]);
        }
        for (int i=0;i<4;i++) for (int j=0;j<4;j++) {
            char nm[32]; snprintf(nm,32,"gram[%d][%d]",i,j);
            pz(nm, gram[i][j]);
        }
        randombytes(r, 8);
        hex("drbg_after_reduce_basis", r, 8);

        // First 4 random samples
        randombytes_init(seed, NULL, 256);
        secret_key_finalize(&sk); secret_key_init(&sk);
        quat_sampling_random_ideal_O0_given_norm(&sk.secret_ideal, &SEC_DEGREE, 1, &QUAT_represent_integer_params, NULL);
        ibz_mat_4x4_finalize(&gram); ibz_mat_4x4_finalize(&red);
        ibz_mat_4x4_init(&gram); ibz_mat_4x4_init(&red);
        quat_lideal_reduce_basis(&red, &gram, &sk.secret_ideal, &QUATALG_PINFTY);
        ibz_t tmp, rmd, an; ibz_init(&tmp); ibz_init(&rmd); ibz_init(&an);
        ibz_mul(&an, &sk.secret_ideal.lattice.denom, &sk.secret_ideal.lattice.denom);
        quat_alg_elem_t na; quat_alg_elem_init(&na);
        for (int iter = 0; iter < 5; iter++) {
            for (int i = 0; i < 4; i++)
                ibz_rand_interval_minm_m(&na.coord[i], QUAT_equiv_bound_coeff);
            quat_qf_eval(&tmp, &gram, &na.coord);
            ibz_div(&tmp, &rmd, &tmp, &an);
            int isp = ibz_probab_prime(&tmp, QUAT_primality_num_iter);
            char nm[64]; snprintf(nm, 64, "iter%d_qf", iter);
            pz(nm, tmp);
            printf("iter%d_isprime = %d\n", iter, isp);
            for (int i=0;i<4;i++){snprintf(nm,64,"iter%d_coord[%d]",iter,i);pz(nm,na.coord[i]);}
            randombytes(r, 8);
            char nh[64]; snprintf(nh,64,"drbg_iter%d",iter);
            hex(nh, r, 8);
        }
    }

    randombytes_init(seed, NULL, 256);
    secret_key_finalize(&sk); secret_key_init(&sk);
    quat_sampling_random_ideal_O0_given_norm(&sk.secret_ideal, &SEC_DEGREE, 1, &QUAT_represent_integer_params, NULL);
    found = found && quat_lideal_prime_norm_reduced_equivalent(
        &sk.secret_ideal, &QUATALG_PINFTY, QUAT_primality_num_iter, QUAT_equiv_bound_coeff);
    printf("after_prime_reduced_found = %d\n", found);
    pz("after_prime_reduced_norm", sk.secret_ideal.norm);

    randombytes(r, 8);
    hex("drbg_after_prime_reduced_next8", r, 8);

    // Re-seed and do REAL protocols_keygen
    randombytes_init(seed, NULL, 256);
    secret_key_finalize(&sk);
    secret_key_init(&sk);
    found = protocols_keygen(&pk, &sk);
    printf("real_keygen_found = %d\n", found);
    pz("real_secret_ideal_norm", sk.secret_ideal.norm);
    randombytes(r, 8);
    hex("drbg_after_real_keygen", r, 8);

    secret_key_finalize(&sk);
    return 0;
}
