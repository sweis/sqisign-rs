fn main() {
    println!("cargo::rustc-check-cfg=cfg(gf5_248_asm)");
    println!("cargo::rustc-check-cfg=cfg(ibz_audit)");
    println!("cargo::rustc-check-cfg=cfg(quaternion_isolation)");

    let gf_fast = std::env::var("CARGO_FEATURE_GF_FAST").is_ok();
    if !gf_fast {
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
            "cargo::warning=gf-fast: target lacks adx+bmi2; using portable backend (set RUSTFLAGS=\"-C target-cpu=native\")"
        );
    }
}
