Experimenting with zkstarks, following https://starkware.co/stark-101/, written in Rust, etc, etc. Mostly finished, just need to finish testing, and implement a random Channel. The stark101 tutorial doesn't cover the proof verifier, so I did my best.

The bulk of the code (that does interesting things) is in prover.rs. It plays out almost 1 to 1 with the stark-101 guide.