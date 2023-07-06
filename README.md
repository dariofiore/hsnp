# Homomorphic Signatures for NP relations


`hsnp` is a Rust library that implements *homomorphic signature schemes for NP relations* (HSNP).

**WARNING:** This is an academic prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Overview
The library currently implements:

1. The HSNP scheme `LHS` from [FT22] for Pedersen commitments to linear functions of signed data.
2. A partial version of the `SPHinx` HSNP scheme for R1CS from [FT22]. The current implementation includes all the building blocks of the generic HSNP construction in [FT22] except that it does not integrate yet `Marlin` zkSNARK.

## Build guide

To compile and execute `hsnp`, you should place `hsnp-poly-commit` (a fork of `ark_poly_commit` available [here](https://github.com/dariofiore/hsnp-poly-commit)) in a sibling directory.

To compile:

```
cargo build
```

To execute the benchmarks, first compile:

```
cargo bench --no-run
```
and then run:
``
./target/release/deps/hsnp_benches-filename HASH_PRECOMP BENCH FIRST LAST SIZE
``
where HASH_PRECOMP= {true|false} and BENCH = {var|hist|mlr|fixed}.
In the case of BENCH=hist SIZE is the size of an histogram bucket.
In the case of BENCH=fixed SIZE is the size of the fixed R1CS.

To execute with memory profiling and to write the timings in a file, run (on Linux):
``
/usr/bin/time -f "\n\nMax memory usage: %M KB." ./target/release/deps/hsnp_benches-filename HASH_PRECOMP BENCH FIRST LAST SIZE > outputfile
``
See [FT22] for more information about the benchmarks.

## License
This code is licensed under either of the following licenses, at your discretion.

- [Apache License Version 2.0](LICENSE-APACHE)
- [MIT License](LICENSE-MIT)

Unless you explicitly state otherwise, any contribution that you submit to this library shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

## Reference paper

[FT22] Dario Fiore, Ida Tucker. Efficient Zero-Knowledge Proofs on Signed Data with Applications to Verifiable Computation on Data Streams. ACM CCS, 2022. https://eprint.iacr.org/2022/1393

## Acknowledgements
This work has received funding by: the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation program under project PICOCRYPT (grant agreement No. 101001283); the Madrid regional government under project BLOQUES (ref. S2018/TCS-4339); a research grant from the Tezos foundation and Nomadic Labs.
