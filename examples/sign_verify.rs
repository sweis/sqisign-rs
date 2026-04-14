//! Generate a SQIsign keypair, sign a message, verify it, and demonstrate
//! that a tampered signature is rejected.
//!
//! Run with: cargo run --release --features lvl1,sign --example sign_verify

use sqisign_rs::common::ctrdrbg;
use sqisign_rs::nistapi::{crypto_sign, crypto_sign_keypair, crypto_sign_open};
use sqisign_rs::params::{CRYPTO_ALGNAME, CRYPTO_BYTES, CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES};

fn main() {
    println!("Algorithm: {CRYPTO_ALGNAME}");
    println!(
        "  pk = {} bytes, sk = {} bytes, sig = {} bytes\n",
        CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES, CRYPTO_BYTES
    );

    // The reference implementation uses an AES-256-CTR-DRBG seeded with 48 bytes.
    // In production this seed would come from the OS RNG.
    let seed = [42u8; 48];
    ctrdrbg::randombytes_init(&seed);

    // Generate a keypair.
    let (pk, sk) = crypto_sign_keypair().expect("keypair generation failed");
    println!("Public key:  {}", hex(&pk));
    println!("Secret key:  {}...\n", hex(&sk[..16]));

    // Sign a message. The output is signature || message.
    let msg = b"Per aspera ad astra";
    let sm = crypto_sign(msg, &sk).expect("signing failed");
    let sig = &sm[..CRYPTO_BYTES];
    println!("Message:     {:?}", std::str::from_utf8(msg).unwrap());
    println!("Signature:   {}\n", hex(sig));

    // Verify the valid signature.
    match crypto_sign_open(&sm, &pk) {
        Ok(recovered) => {
            assert_eq!(recovered, msg);
            println!("Valid signature accepted (recovered message matches).");
        }
        Err(()) => panic!("valid signature was rejected"),
    }

    // Tamper with one byte of the signature and verify it is rejected.
    let mut bad_sm = sm.clone();
    bad_sm[0] ^= 0x01;
    match crypto_sign_open(&bad_sm, &pk) {
        Ok(_) => panic!("tampered signature was accepted"),
        Err(()) => println!("Tampered signature rejected."),
    }

    // Verify against the wrong public key is also rejected.
    let (other_pk, _) = crypto_sign_keypair().expect("keypair generation failed");
    match crypto_sign_open(&sm, &other_pk) {
        Ok(_) => panic!("signature accepted under wrong public key"),
        Err(()) => println!("Signature rejected under unrelated public key."),
    }
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{:02x}", x)).collect()
}
