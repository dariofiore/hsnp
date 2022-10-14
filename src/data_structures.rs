
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
//use ark_ff::Field;

/// Defines the minimal interface for the secret key of any HSNP
/// scheme.
pub trait HSNPSecretKey:
    Clone + core::fmt::Debug // + CanonicalSerialize + CanonicalDeserialize
{

}

/// Defines the minimal interface for public params for any HSNP
/// scheme.
pub trait HSNPUniversalParams:
    Clone + core::fmt::Debug + CanonicalSerialize + CanonicalDeserialize
{
    /// Outputs the maximum degree supported by the committer key.
    fn max_degree(&self) -> usize;
}

/// Defines the minimal interface for labels for any HSNP
/// scheme.
pub trait HSNPLabel:
    Clone + core::fmt::Debug // + CanonicalSerialize + CanonicalDeserialize
{
    
}

/// Defines the minimal interface for signatures for any HSNP
/// scheme.
pub trait HSNPSignature:
    Clone + core::fmt::Debug // + CanonicalSerialize + CanonicalDeserialize
{
    
}

/// Defines the minimal interface for preprocessed evaluation
/// keys for any HSNP scheme.
pub trait PreprocessedEvaluationKey:
    Clone + core::fmt::Debug // + CanonicalSerialize + CanonicalDeserialize
{

}

/// Defines the minimal interface for preprocessed verification
/// keys for any HSNP scheme.
pub trait PreprocessedVerificationKey:
    Clone + core::fmt::Debug // + CanonicalSerialize + CanonicalDeserialize
{

}

///defines the interface for statements of NP relations
pub trait RelationStatement:
    Clone + core::fmt::Debug
{

}


///defines the interface for signed inputs of NP relations
pub trait RelationSignedInput:
    Clone + core::fmt::Debug
{

}

///defines the interface for witnesses of NP relations
pub trait RelationWitness:
    Clone + core::fmt::Debug
{

}

/// Defines the minimal interface for the relations to be proven in
/// the HSNP scheme.
pub trait Relation:
    Clone + core::fmt::Debug
{
    ///The type of statements for this relation
    type Statement: RelationStatement;
    ///the type of signed inputs for this relation
    type SignedInput: RelationSignedInput;
    ///the type of witnesses for this relation
    type Witness: RelationWitness;

    ///checks whether the relation holds
    fn check(
        st: &Self::Statement,
        sig_inp: &Self::SignedInput,
        wit: &Self::Witness,
    ) -> bool;
}

/// Defines the minimal interface for the proofs generated in
/// the HSNP scheme.
pub trait HSNPProof:
    Clone + core::fmt::Debug
{

}