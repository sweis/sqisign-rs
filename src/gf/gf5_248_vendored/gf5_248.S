.intel_syntax noprefix

#ifdef __APPLE__
.section __TEXT,__const
#else
.section .rodata
#endif

p:        .quad 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x04ffffffffffffff
p2:       .quad 0xfffffffffffffffe, 0xffffffffffffffff, 0xffffffffffffffff, 0x09ffffffffffffff
p_plus_1: .quad 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0500000000000000

#if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",@progbits
#endif

.text
.p2align 4,,15

#ifdef __APPLE__
.global _fp_add_asm
_fp_add_asm:
#else
.global fp_add_asm
fp_add_asm:
#endif
  xor    rax, rax
  mov    r8, [rsi]
  mov    r9, [rsi+8]
  mov    r10, [rsi+16]
  mov    r11, [rsi+24]
  add    r8, [rdx]
  adc    r9, [rdx+8]
  adc    r10, [rdx+16]
  adc    r11, [rdx+24]
  mov    rax, r11
  shr    rax, 59
  neg    rax
  mov    rdx, [rip+p+24]
  and    rdx, rax
  sub    r8, rax
  sbb    r9, rax
  sbb    r10, rax
  sbb    r11, rdx

  mov    rax, r11
  shr    rax, 59
  neg    rax
  mov    rdx, [rip+p+24]
  and    rdx, rax
  sub    r8, rax
  sbb    r9, rax
  sbb    r10, rax
  sbb    r11, rdx
  mov    [rdi], r8
  mov    [rdi+8], r9
  mov    [rdi+16], r10
  mov    [rdi+24], r11
  ret

#ifdef __APPLE__
.global _fp_sub_asm
_fp_sub_asm:
#else
.global fp_sub_asm
fp_sub_asm:
#endif
  xor    rax, rax
  mov    r8, [rsi]
  mov    r9, [rsi+8]
  mov    r10, [rsi+16]
  mov    r11, [rsi+24]
  sub    r8, [rdx]
  sbb    r9, [rdx+8]
  sbb    r10, [rdx+16]
  sbb    r11, [rdx+24]
  sbb    rax, 0

  mov    rdx, [rip+p+24]
  and    rdx, rax
  add    r8, rax
  adc    r9, rax
  adc    r10, rax
  adc    r11, rdx

  mov    rax, r11
  sar    rax, 59
  mov    rdx, [rip+p+24]
  and    rdx, rax
  add    r8, rax
  adc    r9, rax
  adc    r10, rax
  adc    r11, rdx

  mov    [rdi], r8
  mov    [rdi+8], r9
  mov    [rdi+16], r10
  mov    [rdi+24], r11
  ret

///////////////////////////////////////////////////////////////// MACROS
// z = a x bi + z
// Inputs: base memory pointer M1 (a),
//         bi pre-stored in rdx,
//         accumulator z in [Z0:Z4]
// Output: [Z0:Z4]
// Temps:  regs T0:T1
/////////////////////////////////////////////////////////////////
.macro MULADD64x256 M1, Z0, Z1, Z2, Z3, Z4, T0, T1, C
    mulx   \T0, \T1, \M1     // A0*B0
    xor    \C, \C
    adox   \Z0, \T1
    adox   \Z1, \T0
    mulx   \T0, \T1, 8\M1    // A0*B1
    adcx   \Z1, \T1
    adox   \Z2, \T0
    mulx   \T0, \T1, 16\M1   // A0*B2
    adcx   \Z2, \T1
    adox   \Z3, \T0
    mulx   \T0, \T1, 24\M1   // A0*B3
    adcx   \Z3, \T1
    adox   \Z4, \T0
    adc    \Z4, 0
.endm

.macro MULADD64x64 M1, Z0, Z1, Z2, Z3, T0, T1
    mulx   \T0, \T1, \M1     // A0*B0
    xor    rax, rax
    adox   \Z2, \T1
    adox   \Z3, \T0
.endm

///////////////////////////////////////////////////////////////// MACRO
// z = a x b (mod p)
// Inputs: base memory pointers M0 (a), M1 (b)
//         bi pre-stored in rdx,
//         accumulator z in [Z0:Z4], pre-stores a0 x b
// Output: [Z0:Z4]
// Temps:  regs T0:T1
/////////////////////////////////////////////////////////////////
.macro FPMUL256x256 M0, M1, Z0, Z1, Z2, Z3, Z4, T0, T1
    // [Z1:Z4] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, \Z0                 // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], \Z1, \Z2, \Z3, \Z4, \T0, \T1

    // [Z1:Z4, Z0] <- z = a01 x a1 + z
    mov    rdx, 8\M0
    MULADD64x256 \M1, \Z1, \Z2, \Z3, \Z4, \Z0, \T0, \T1, \Z0
    // [Z2:Z4, Z0] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, \Z1                 // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], \Z2, \Z3, \Z4, \Z0, \T0, \T1

    // [Z2:Z4, Z0:Z1] <- z = a02 x a1 + z
    mov    rdx, 16\M0
    MULADD64x256 \M1, \Z2, \Z3, \Z4, \Z0, \Z1, \T0, \T1, \Z1
    // [Z3:Z4, Z0:Z1] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, \Z2                // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], \Z3, \Z4, \Z0, \Z1, \T0, \T1

    // [Z3:Z4, Z0:Z2] <- z = a03 x a1 + z
    mov    rdx, 24\M0
    MULADD64x256 \M1, \Z3, \Z4, \Z0, \Z1, \Z2, \T0, \T1, \Z2
    // [Z4, Z0:Z2] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, \Z3                // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], \Z4, \Z0, \Z1, \Z2, \T0, \T1
.endm

#ifdef __APPLE__
.global _fp_mul_asm
_fp_mul_asm:
#else
.global fp_mul_asm
fp_mul_asm:
#endif
    push   r12
    push   r13
    push   r14
    mov    rcx, rdx

    // [r8:r12] <- z = a x b0
    mov    rdx, [rcx]
    mulx   r9, r8, [rsi]
    xor    rax, rax
    mulx   r10, r11, [rsi+8]
    adox   r9, r11
    mulx   r11, r12, [rsi+16]
    adox   r10, r12
    mulx   r12, r13, [rsi+24]
    adox   r11, r13
    adox   r12, rax

	FPMUL256x256 [rcx], [rsi], r8, r9, r10, r11, r12, r13, r14

    mov    [rdi], r12
    mov    [rdi+8], r8
    mov    [rdi+16], r9
    mov    [rdi+24], r10
    pop    r14
    pop    r13
    pop    r12
    ret

//***********************************************************************
//  Squaring in GF(p)
//  Operation: c [rdi] = a^2 (mod p), a in [rsi]
//
//  Schoolbook square (10 MULX vs 16 for fp_mul) followed by the same
//  special-form Montgomery reduction (4 MULX) as fp_mul. Doubling and
//  diagonal addition are interleaved on independent CF/OF chains
//  (ADCX for doubling, ADOX for diagonals).
//
//  Output is bit-identical to fp_mul_asm(c, a, a) for any input < 2^252.
//***********************************************************************
#ifdef __APPLE__
.global _fp_sqr_asm
_fp_sqr_asm:
#else
.global fp_sqr_asm
fp_sqr_asm:
#endif
    push   r12
    push   r13
    push   r14
    push   r15
    push   rbx

    // --- Phase 1a: off-diagonal products into e1..e6 = r8..r13 ---
    mov    rdx, [rsi]
    mulx   r9, r8, [rsi+8]       // a0*a1 -> e1:e2
    mulx   r11, r10, [rsi+24]    // a0*a3 -> e3:e4
    mov    rdx, [rsi+16]
    mulx   r13, r12, [rsi+24]    // a2*a3 -> e5:e6

    mov    rdx, [rsi]
    mulx   rcx, rbx, [rsi+16]    // a0*a2
    add    r9,  rbx
    adc    r10, rcx
    mov    rdx, [rsi+8]
    mulx   rcx, rbx, [rsi+24]    // a1*a3
    adc    r11, rbx
    adc    r12, rcx
    adc    r13, 0
    mulx   rcx, rbx, [rsi+16]    // a1*a2 (rdx still a1)
    add    r10, rbx
    adc    r11, rcx
    adc    r12, 0
    adc    r13, 0

    // --- Phase 1b: 2*(off-diag) + diag, via dual carry chains ---
    // ADCX (CF) doubles e1..e6 in place; ADOX (OF) accumulates diagonals.
    xor    eax, eax              // CF=OF=0
    mov    rdx, [rsi]
    mulx   rbx, r15, rdx         // a0^2 -> e0=r15 : hi=rbx
    adcx   r8, r8
    adox   r8, rbx
    mov    rdx, [rsi+8]
    mulx   rcx, rbx, rdx         // a1^2
    adcx   r9, r9
    adox   r9, rbx
    adcx   r10, r10
    adox   r10, rcx
    mov    rdx, [rsi+16]
    mulx   rcx, rbx, rdx         // a2^2
    adcx   r11, r11
    adox   r11, rbx
    adcx   r12, r12
    adox   r12, rcx
    mov    rdx, [rsi+24]
    mulx   rcx, rbx, rdx         // a3^2
    adcx   r13, r13
    adox   r13, rbx
    mov    r14, 0                // mov: flags untouched
    adcx   r14, r14              // r14 = CF (top bit of doubled e6); CF cleared
    adox   r14, rcx              // e7 = top + hi(a3^2) + OF

    // e[0..7] = r15, r8, r9, r10, r11, r12, r13, r14

    // --- Phase 2: 4-step Montgomery reduction ---
    // For each m = e_i, add m*5<<56 into e[i+3], e[i+4], propagate carry to e7.
    mov    rdx, r15
    mulx   rcx, rbx, [rip+p_plus_1+24]
    add    r10, rbx
    adc    r11, rcx
    adc    r12, 0
    adc    r13, 0
    adc    r14, 0

    mov    rdx, r8
    mulx   rcx, rbx, [rip+p_plus_1+24]
    add    r11, rbx
    adc    r12, rcx
    adc    r13, 0
    adc    r14, 0

    mov    rdx, r9
    mulx   rcx, rbx, [rip+p_plus_1+24]
    add    r12, rbx
    adc    r13, rcx
    adc    r14, 0

    mov    rdx, r10              // e3 was updated by step 0
    mulx   rcx, rbx, [rip+p_plus_1+24]
    add    r13, rbx
    adc    r14, rcx

    mov    [rdi],    r11
    mov    [rdi+8],  r12
    mov    [rdi+16], r13
    mov    [rdi+24], r14

    pop    rbx
    pop    r15
    pop    r14
    pop    r13
    pop    r12
    ret

//***********************************************************************
//  Multiplication in GF(p^2), non-complex part
//  Operation: c [rdi] = a0 x b0 - a1 x b1
//  Inputs: a = [a1, a0] stored in [rsi]
//          b = [b1, b0] stored in [rdx]
//  Output: c stored in [rdi]
//***********************************************************************
#ifdef __APPLE__
.global _fp_dop_asm
_fp_dop_asm:
#else
.global fp_dop_asm
fp_dop_asm:
#endif
    push   r12
    push   r13
    push   r14
    mov    rcx, rdx

	// [rdi0:3] <- 2p - b1
	mov    r8, [rip+p2]
	mov    r9, [rip+p2+8]
	mov    r10, r9
	mov    r11, [rip+p2+24]
	mov    rax, [rcx+32]
	mov    rdx, [rcx+40]
	sub    r8, rax
	sbb    r9, rdx
	mov    rax, [rcx+48]
	mov    rdx, [rcx+56]
	sbb    r10, rax
	sbb    r11, rdx
	mov    [rdi], r8
	mov    [rdi+8], r9
	mov    [rdi+16], r10
	mov    [rdi+24], r11

    // [r8:r12] <- z = a0 x b00 - a1 x b10
    mov    rdx, [rcx]
    mulx   r9, r8, [rsi]
    xor    rax, rax
    mulx   r10, r11, [rsi+8]
    adox   r9, r11
    mulx   r11, r12, [rsi+16]
    adox   r10, r12
    mulx   r12, r13, [rsi+24]
    adox   r11, r13
    adox   r12, rax

    mov    rdx, [rdi]
    MULADD64x256 [rsi+32], r8, r9, r10, r11, r12, r13, r14, rax
    // [r9:r12] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, r8                 // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], r9, r10, r11, r12, r13, r14

    // [r9:r12, r8] <- z = a0 x b01 - a1 x b11 + z
    mov    rdx, [rcx+8]
    MULADD64x256 [rsi], r9, r10, r11, r12, r8, r13, r14, r8
    mov    rdx, [rdi+8]
    MULADD64x256 [rsi+32], r9, r10, r11, r12, r8, r13, r14, rax
    // [r10:r12, r8] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, r9                 // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], r10, r11, r12, r8, r13, r14

    // [r10:r12, r8:r9] <- z = a0 x b02 - a1 x b12 + z
    mov    rdx, [rcx+16]
    MULADD64x256 [rsi], r10, r11, r12, r8, r9, r13, r14, r9
    mov    rdx, [rdi+16]
    MULADD64x256 [rsi+32], r10, r11, r12, r8, r9, r13, r14, rax
    // [r11:r12, r8:r9] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, r10                // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], r11, r12, r8, r9, r13, r14

    // [r11:r12, r8:r10] <- z = a0 x b03 - a1 x b13 + z
    mov    rdx, [rcx+24]
    MULADD64x256 [rsi], r11, r12, r8, r9, r10, r13, r14, r10
    mov    rdx, [rdi+24]
    MULADD64x256 [rsi+32], r11, r12, r8, r9, r10, r13, r14, rax
    // [r12, r8:r10] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, r11                // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], r12, r8, r9, r10, r13, r14

    mov    [rdi], r12
    mov    [rdi+8], r8
    mov    [rdi+16], r9
    mov    [rdi+24], r10
    pop    r14
    pop    r13
    pop    r12
    ret

//***********************************************************************
//  Multiplication in GF(p^2), complex part
//  Operation: c [rdi] = a0 x b1 + a1 x b0
//  Inputs: a = [a1, a0] stored in [rsi]
//          b = [b1, b0] stored in [rdx]
//  Output: c stored in [rdi]
//***********************************************************************
#ifdef __APPLE__
.global _fp_sop_asm
_fp_sop_asm:
#else
.global fp_sop_asm
fp_sop_asm:
#endif
    push   r12
    push   r13
    push   r14
    mov    rcx, rdx

    // [r8:r12] <- z = a0 x b10 + a1 x b00
    mov    rdx, [rcx+32]
    mulx   r9, r8, [rsi]
    xor    rax, rax
    mulx   r10, r11, [rsi+8]
    adox   r9, r11
    mulx   r11, r12, [rsi+16]
    adox   r10, r12
    mulx   r12, r13, [rsi+24]
    adox   r11, r13
    adox   r12, rax

    mov    rdx, [rcx]
    MULADD64x256 [rsi+32], r8, r9, r10, r11, r12, r13, r14, rax
    // [r9:r12] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, r8                 // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], r9, r10, r11, r12, r13, r14

    // [r9:r12, r8] <- z = a0 x b01 - a1 x b11 + z
    mov    rdx, [rcx+40]
    MULADD64x256 [rsi], r9, r10, r11, r12, r8, r13, r14, r8
    mov    rdx, [rcx+8]
    MULADD64x256 [rsi+32], r9, r10, r11, r12, r8, r13, r14, rax
    // [r10:r12, r8] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, r9                 // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], r10, r11, r12, r8, r13, r14

    // [r10:r12, r8:r9] <- z = a0 x b02 - a1 x b12 + z
    mov    rdx, [rcx+48]
    MULADD64x256 [rsi], r10, r11, r12, r8, r9, r13, r14, r9
    mov    rdx, [rcx+16]
    MULADD64x256 [rsi+32], r10, r11, r12, r8, r9, r13, r14, rax
    // [r11:r12, r8:r9] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, r10                // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], r11, r12, r8, r9, r13, r14

    // [r11:r12, r8:r10] <- z = a0 x b03 - a1 x b13 + z
    mov    rdx, [rcx+56]
    MULADD64x256 [rsi], r11, r12, r8, r9, r10, r13, r14, r10
    mov    rdx, [rcx+24]
    MULADD64x256 [rsi+32], r11, r12, r8, r9, r10, r13, r14, rax
    // [r12, r8:r10] <- z = (z0 x p_plus_1 + z)/2^64
    mov    rdx, r11                // rdx <- z0
    MULADD64x64 [rip+p_plus_1+24], r12, r8, r9, r10, r13, r14

    mov    [rdi], r12
    mov    [rdi+8], r8
    mov    [rdi+16], r9
    mov    [rdi+24], r10
    pop    r14
    pop    r13
    pop    r12
    ret

//***********************************************************************
//  Squaring in GF(p^2), non-complex part
//  Operation: c [rdi] = (a0+a1) x (a0-a1)
//  Inputs: a = [a1, a0] stored in [rsi]
//  Output: c stored in [rdi]
//***********************************************************************
#ifdef __APPLE__
.global _fp2_sqr_re_asm
_fp2_sqr_re_asm:
#else
.global fp2_sqr_re_asm
fp2_sqr_re_asm:
#endif
    push   r12
    push   r13

	// a0 + a1
	mov    rdx, [rsi]
	mov    r9, [rsi+8]
	mov    r10, [rsi+16]
	mov    r11, [rsi+24]
	add    rdx, [rsi+32]
	adc    r9, [rsi+40]
	adc    r10, [rsi+48]
	adc    r11, [rsi+56]
	mov    [rdi], rdx
	mov    [rdi+8], r9
	mov    [rdi+16], r10
	mov    [rdi+24], r11

	// a0 - a1 + 2p
	mov    r8, [rsi]
	mov    r10, [rsi+8]
	mov    r12, [rsi+16]
	mov    r13, [rsi+24]
	sub    r8, [rsi+32]
	sbb    r10, [rsi+40]
	sbb    r12, [rsi+48]
	sbb    r13, [rsi+56]
	mov    rax, [rip+p2]
	add    r8, rax
	mov    rax, [rip+p2+8]
	adc    r10, rax
	adc    r12, rax
	adc    r13, [rip+p2+24]
	mov    [rdi+32], r8
	mov    [rdi+40], r10
	mov    [rdi+48], r12
	mov    [rdi+56], r13

    // [r8:r12] <- z = a00 x a1
    mulx   r9, r8, r8
    xor    rax, rax
    mulx   r10, r11, r10
    adox   r9, r11
    mulx   r11, r12, r12
    adox   r10, r12
    mulx   r12, r13, r13
    adox   r11, r13
    adox   r12, rax

    FPMUL256x256 [rdi], [rdi+32], r8, r9, r10, r11, r12, r13, rcx

    mov    [rdi], r12
    mov    [rdi+8], r8
    mov    [rdi+16], r9
    mov    [rdi+24], r10
    pop    r13
    pop    r12
    ret

//***********************************************************************
//  Squaring in GF(p^2), complex part
//  Operation: c [rdi] = 2a0 x a1
//  Inputs: a = [a1, a0] stored in [reg_p1]
//  Output: c stored in [rdi]
//***********************************************************************
#ifdef __APPLE__
.global _fp2_sqr_im_asm
_fp2_sqr_im_asm:
#else
.global fp2_sqr_im_asm
fp2_sqr_im_asm:
#endif
    push   r12
    push   r13

	mov    rdx, [rsi]
	mov    r9, [rsi+8]
	mov    r10, [rsi+16]
	mov    r11, [rsi+24]
	add    rdx, rdx
	adc    r9, r9
	adc    r10, r10
	adc    r11, r11
	sub    rsp, 32
	mov    [rsp+8], r9
	mov    [rsp+16], r10
	mov    [rsp+24], r11

    // [r8:r12] <- z = a00 x a1
    mulx   r9, r8, [rsi+32]
    xor    rax, rax
    mulx   r10, r11, [rsi+40]
    adox   r9, r11
    mulx   r11, r12, [rsi+48]
    adox   r10, r12
    mulx   r12, r13, [rsi+56]
    adox   r11, r13
    adox   r12, rax

	FPMUL256x256 [rsp], [rsi+32], r8, r9, r10, r11, r12, r13, rcx
	add    rsp, 32

    mov    [rdi], r12
    mov    [rdi+8], r8
    mov    [rdi+16], r9
    mov    [rdi+24], r10
    pop    r13
    pop    r12
    ret
