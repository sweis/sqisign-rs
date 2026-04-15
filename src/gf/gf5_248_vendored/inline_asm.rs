//! Inline `asm!()` versions of `fp_mul`/`fp_sqr` for the GF5_248 prime
//! (5·2²⁴⁸ − 1). Algorithmically identical to `gf5_248.S` but inlineable,
//! letting LLVM schedule across field operations.
//!
//! Requires x86_64 with BMI2 (`mulx`) and ADX (`adcx`/`adox`).
#![allow(unsafe_code)]

use core::arch::asm;

const P_PLUS_1_HI: u64 = 0x0500_0000_0000_0000; // (p+1) / 2¹⁹²

/// `r ← a · b · R⁻¹ mod p` (Montgomery). Interleaved CIOS, 20 MULX.
#[inline(always)]
pub fn fp_mul(r: &mut [u64; 4], a: &[u64; 4], b: &[u64; 4]) {
    let (o0, o1, o2, o3): (u64, u64, u64, u64);
    unsafe {
        asm!(
            // Row 0: z = a × b₀  → [r8:r12]
            "mov  rdx, [{b}]",
            "mulx r9, r8, [{a}]",
            "xor  rax, rax",
            "mulx r10, r11, [{a}+8]",
            "adox r9, r11",
            "mulx r11, r12, [{a}+16]",
            "adox r10, r12",
            "mulx r12, r13, [{a}+24]",
            "adox r11, r13",
            "adox r12, rax",

            // ---- reduce₀ ----
            "mov  rdx, r8",
            "mulx r13, r14, {pp1}",
            "xor  rax, rax",
            "adox r11, r14",
            "adox r12, r13",
            // ---- row 1 into [r9:r12,r8] ----
            "mov  rdx, [{b}+8]",
            "mulx r13, r14, [{a}]",
            "xor  r8, r8",
            "adox r9, r14",
            "adox r10, r13",
            "mulx r13, r14, [{a}+8]",
            "adcx r10, r14",
            "adox r11, r13",
            "mulx r13, r14, [{a}+16]",
            "adcx r11, r14",
            "adox r12, r13",
            "mulx r13, r14, [{a}+24]",
            "adcx r12, r14",
            "adox r8, r13",
            "adc  r8, 0",
            // ---- reduce₁ ----
            "mov  rdx, r9",
            "mulx r13, r14, {pp1}",
            "xor  rax, rax",
            "adox r12, r14",
            "adox r8, r13",
            // ---- row 2 into [r10:r12,r8,r9] ----
            "mov  rdx, [{b}+16]",
            "mulx r13, r14, [{a}]",
            "xor  r9, r9",
            "adox r10, r14",
            "adox r11, r13",
            "mulx r13, r14, [{a}+8]",
            "adcx r11, r14",
            "adox r12, r13",
            "mulx r13, r14, [{a}+16]",
            "adcx r12, r14",
            "adox r8, r13",
            "mulx r13, r14, [{a}+24]",
            "adcx r8, r14",
            "adox r9, r13",
            "adc  r9, 0",
            // ---- reduce₂ ----
            "mov  rdx, r10",
            "mulx r13, r14, {pp1}",
            "xor  rax, rax",
            "adox r8, r14",
            "adox r9, r13",
            // ---- row 3 into [r11:r12,r8,r9,r10] ----
            "mov  rdx, [{b}+24]",
            "mulx r13, r14, [{a}]",
            "xor  r10, r10",
            "adox r11, r14",
            "adox r12, r13",
            "mulx r13, r14, [{a}+8]",
            "adcx r12, r14",
            "adox r8, r13",
            "mulx r13, r14, [{a}+16]",
            "adcx r8, r14",
            "adox r9, r13",
            "mulx r13, r14, [{a}+24]",
            "adcx r9, r14",
            "adox r10, r13",
            "adc  r10, 0",
            // ---- reduce₃ ----
            "mov  rdx, r11",
            "mulx r13, r14, {pp1}",
            "xor  rax, rax",
            "adox r9, r14",
            "adox r10, r13",

            a = in(reg) a.as_ptr(),
            b = in(reg) b.as_ptr(),
            pp1 = in(reg) P_PLUS_1_HI,
            out("rax") _, out("rdx") _,
            out("r8") o1, out("r9") o2, out("r10") o3, out("r11") _,
            out("r12") o0, out("r13") _, out("r14") _,
            options(nostack, pure, readonly),
        );
    }
    r[0] = o0;
    r[1] = o1;
    r[2] = o2;
    r[3] = o3;
}

/// `r ← a² · R⁻¹ mod p` (Montgomery). 14 MULX (10 square + 4 reduce).
/// Bit-identical to `fp_mul(r, a, a)` for inputs < 2²⁵².
#[inline(always)]
pub fn fp_sqr(r: &mut [u64; 4], a: &[u64; 4]) {
    let (o0, o1, o2, o3): (u64, u64, u64, u64);
    unsafe {
        asm!(
            // ---- Phase 1a: off-diagonal products into e1..e6 = r8..r13 ----
            "mov  rdx, [{a}]",
            "mulx r9, r8, [{a}+8]",       // a0*a1 -> e1:e2
            "mulx r11, r10, [{a}+24]",    // a0*a3 -> e3:e4
            "mov  rdx, [{a}+16]",
            "mulx r13, r12, [{a}+24]",    // a2*a3 -> e5:e6

            "mov  rdx, [{a}]",
            "mulx rcx, rax, [{a}+16]",    // a0*a2
            "add  r9,  rax",
            "adc  r10, rcx",
            "mov  rdx, [{a}+8]",
            "mulx rcx, rax, [{a}+24]",    // a1*a3
            "adc  r11, rax",
            "adc  r12, rcx",
            "adc  r13, 0",
            "mulx rcx, rax, [{a}+16]",    // a1*a2 (rdx still a1)
            "add  r10, rax",
            "adc  r11, rcx",
            "adc  r12, 0",
            "adc  r13, 0",

            // ---- Phase 1b: 2*(off-diag) + diag, ADCX doubles || ADOX diag ----
            "xor  r15d, r15d",            // clear CF/OF; r15 overwritten next
            "mov  rdx, [{a}]",
            "mulx rax, r15, rdx",         // a0² -> e0=r15 : hi=rax
            "adcx r8, r8",
            "adox r8, rax",
            "mov  rdx, [{a}+8]",
            "mulx rcx, rax, rdx",         // a1²
            "adcx r9, r9",
            "adox r9, rax",
            "adcx r10, r10",
            "adox r10, rcx",
            "mov  rdx, [{a}+16]",
            "mulx rcx, rax, rdx",         // a2²
            "adcx r11, r11",
            "adox r11, rax",
            "adcx r12, r12",
            "adox r12, rcx",
            "mov  rdx, [{a}+24]",
            "mulx rcx, rax, rdx",         // a3²
            "adcx r13, r13",
            "adox r13, rax",
            "mov  r14, 0",
            "adcx r14, r14",
            "adox r14, rcx",

            // e[0..7] = r15, r8, r9, r10, r11, r12, r13, r14

            // ---- Phase 2: 4-step Montgomery reduction ----
            // Load (p+1)>>192 into r15 once e0 has been moved to rdx.
            "mov  rdx, r15",
            "mov  r15, {pp1}",
            "mulx rcx, rax, r15",
            "add  r10, rax",
            "adc  r11, rcx",
            "adc  r12, 0",
            "adc  r13, 0",
            "adc  r14, 0",

            "mov  rdx, r8",
            "mulx rcx, rax, r15",
            "add  r11, rax",
            "adc  r12, rcx",
            "adc  r13, 0",
            "adc  r14, 0",

            "mov  rdx, r9",
            "mulx rcx, rax, r15",
            "add  r12, rax",
            "adc  r13, rcx",
            "adc  r14, 0",

            "mov  rdx, r10",
            "mulx rcx, rax, r15",
            "add  r13, rax",
            "adc  r14, rcx",

            a = in(reg) a.as_ptr(),
            pp1 = const P_PLUS_1_HI,
            out("rax") _, out("rcx") _, out("rdx") _,
            out("r8") _, out("r9") _, out("r10") _,
            out("r11") o0, out("r12") o1, out("r13") o2, out("r14") o3,
            out("r15") _,
            options(nostack, pure, readonly),
        );
    }
    r[0] = o0;
    r[1] = o1;
    r[2] = o2;
    r[3] = o3;
}
