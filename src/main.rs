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
    // Abstracts the interactive verifier
    let channel = Channel::new();

    // Generates a proof, using the channel to provide data
    let proof = generate_proof(channel);

    // Verify the proof (this will panic if anything goes wrong, proper error handling comes later)
    proof.verify();

    // Yay, we did it. Print proof size.
    println!("Proof succes! Proof size: {:?}", proof.size());
}
