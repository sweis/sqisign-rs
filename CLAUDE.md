This project is going to be a Rust implementation of SQISign.

For general information on SQISign in general: https://sqisign.org/

The specification is here: https://sqisign.org/spec/sqisign-20250707.pdf

The source in C is here: https://github.com/SQIsign/the-sqisign

Another reference implementation in Sage is here: https://github.com/KULeuven-COSIC/sqisign_prism_v2

The plan should be:
1. Import any test vectors from C and the specification and use those.
2. Port the C code to Rust
3. Cross-check the specification and validate that the C and Rust implementations are implementing it correctly. If there is something that looks incorrect, fix it in Rust and document what is wrong with C.
4. If there are additional test vectors that do not exist and would help confirm correctness, add them. 
5. Once we have parity with Rust and C, do a cross-compatibility check with Sage. Do not take the Sage implementation as authoriatative and note any discrepancies.
6. Once we are confident that the Rust implementation is correct, look for security and performance improvements that will make it superior to the C implementation.
