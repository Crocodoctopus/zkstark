mod channel;
mod field;
mod merkle;
mod polynomial;
mod proof;
mod prover;

use channel::Channel;
use prover::generate_proof;

// Represents an element of a prime field
// All math is done mod 3221225473
type F = field::Gf<3221225473>;

fn main() {
    use std::time::Instant;

    // Abstracts the interactive verifier
    let channel = Channel::new();

    // Generates a proof, using the channel to provide data
    let start = Instant::now();
    let proof = generate_proof(channel);
    println!("Prover runtime: {:?}", Instant::now().duration_since(start));

    // Verify the proof (this will panic if anything goes wrong, proper error handling comes later)
    let start = Instant::now();
    proof.verify();
    println!(
        "Verifier runtime: {:?}",
        Instant::now().duration_since(start)
    );

    // Yay, we did it. Print proof size.
    println!("Proof size: {:?}", proof.size());
}
