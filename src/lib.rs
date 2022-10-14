#![cfg_attr(not(feature = "std"), no_std)]
//! A crate for HSNP schemes.
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, private_in_public, variant_size_differences)]
#![deny(stable_features, unreachable_pub, non_shorthand_field_patterns)]
#![deny(unused_attributes, unused_mut)]
#![deny(missing_docs)]
#![deny(unused_imports)]
#![deny(renamed_and_removed_lints, stable_features, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, unused_must_use, const_err)]
#![forbid(unsafe_code)]

#[macro_use]
extern crate ark_std;
#[macro_use]
extern crate derivative;

use ark_ff::PrimeField;
use ark_std::rand::RngCore;
//use ark_poly_commit::Error;

/// Data structures used by an HSNP scheme.
pub mod data_structures;
pub use data_structures::*;

/// Errors pertaining to query sets.
pub mod error;
pub use error::*;


///Implements the FT22 HSNP scheme for R1CS relations
pub mod sphinx_hsnp;

///Implements the FT22 HSNP scheme for relations verifying
///  Pedersen-committed outputs of linear functions
pub mod lhs_ped;

///Describes the interface of an HSNP scheme for relations over
///signed elements of a prime finite field.
pub trait HSNP<F: PrimeField, Rel: Relation>:
    Sized 
{
    ///The secret key of the HSNP scheme
    type SecretKey: HSNPSecretKey;

    ///The Universal parameters of the HSNP scheme
    type UniversalParams: HSNPUniversalParams;

    ///The preprocessed evaluation key for a specific relation
    type EvKeyR: PreprocessedEvaluationKey;

    ///The preprocessed verification key for a specific relation
    type VkR: PreprocessedVerificationKey;

    ///The signature of the HSNP scheme
    type Signature: HSNPSignature;

    ///The label of the HSNP scheme
    type Label: HSNPLabel;

    
    ///The proof produced by the HSNP scheme
    type Proof: HSNPProof;

    /// The error type for the scheme.
    //type Error: ark_std::error::Error + From<Error>;
    type Error: ark_std::error::Error;



    /// Constructs public parameters  when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme. `max_inputs` specifies the maximum number
    /// of signed inputs
    fn setup<R: RngCore>(
        max_degree: usize,
        max_inputs: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error>;

    ///Specializes the universal parameters to a given relation
    fn preprocess(
        pp: &Self::UniversalParams,
        rel: &Rel,
    ) -> (Self::EvKeyR, Self::VkR);

    ///signing algorithm
    fn sign<R: RngCore>(
        sk: &Self::SecretKey,
        lab: &Self::Label,
        msg: &F,
    ) -> Self::Signature;

    ///HSNP evaluation algorithm
    fn eval<R: RngCore>(
        evk: &Self::EvKeyR,
        rel: &Rel,
        stmt: &Rel::Statement,
        sigs: &[Self::Signature],
    ) -> Self::Proof;

    ///HSNP verification algorithm
    fn ver(
        vk: &Self::VkR,
        labels: &[Self::Label],
        stmt: &Rel::Statement,
        proof: &[Self::Proof],
    ) -> Result<bool, Self::Error>;

}