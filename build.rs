fn main() {
    println!("cargo::rustc-check-cfg=cfg(gf5_248_asm)");
    println!("cargo::rustc-check-cfg=cfg(ibz_audit)");
    println!("cargo::rustc-check-cfg=cfg(quaternion_isolation)");

    // The vendored asm GF backend is the default at lvl1. lvl3/lvl5 and the
    // gf-portable opt-out use modarith (no asm).
    let lvl1 = std::env::var("CARGO_FEATURE_LVL1").is_ok();
    let lvl3 = std::env::var("CARGO_FEATURE_LVL3").is_ok();
    let lvl5 = std::env::var("CARGO_FEATURE_LVL5").is_ok();
    let portable = std::env::var("CARGO_FEATURE_GF_PORTABLE").is_ok();
    if !lvl1 || lvl3 || lvl5 || portable {
        return;
    }

    let arch = std::env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let features = std::env::var("CARGO_CFG_TARGET_FEATURE").unwrap_or_default();
    let has_adx_bmi2 = arch == "x86_64" && features.contains("adx") && features.contains("bmi2");

    if has_adx_bmi2 {
        cc::Build::new()
            .file("src/gf/gf5_248_vendored/gf5_248.S")
            .flag("-O3")
            .compile("gf5_248");
        println!("cargo::rustc-cfg=gf5_248_asm");
    } else {
        println!(
            "cargo::warning=lvl1 GF: target lacks adx+bmi2; using portable fallback (set RUSTFLAGS=\"-C target-cpu=native\" for asm)"
        );
    }
}
