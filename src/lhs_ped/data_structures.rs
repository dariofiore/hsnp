use ark_poly_commit::hsnp_pc;
use ark_ec::PairingEngine;
use ark_ff::ToBytes;
use ark_std::io::Write;

/// Public statement of the relation
pub type PedLinRelationStatement<E> = hsnp_pc::Commitment<E>;

/// Witness of the relation
pub type PedLinRelationWitness<E> = <E as PairingEngine>::Fr;

/// Signed input of the relation
pub type PedLinRelationSignedInput<E> = <E as PairingEngine>::Fr;

///Description of the universal relation
pub struct UniversalPedLinRelation<E: PairingEngine> {
    ///(short) commitment key
    pub ck: hsnp_pc::ScalarComKey<E>,
        
}

///Description of the relation R_f(c, z, r)
pub struct PedLinRelation<E: PairingEngine> {
    
    ///length of the signed input
    pub signed_input_size: usize,
    //pub stmt_size: usize,

    ///the linear function expresses as a vector of coefficients
    pub f: Vec<E::Fr>,
    
}

///Implements the public key of the HSNP scheme for committed linear functions from [FT22]
#[derive(Derivative)]
#[derivative(
    Clone(bound = ""),
)]
pub struct HSNPLinPublicKey<E: PairingEngine> {
    ///the commitment key for the HSNP_PC KZG10-based commitment to scalars
    pub ck: hsnp_pc::ScalarComKey<E>,
    ///element of the public key
    pub gamma1: E::G1Affine,
    ///element of the public key
    pub gamma2: E::G2Affine,
    ///element of the public key
    pub b_g: E::G1Affine,
}

///Implements the secret key of the HSNP scheme for committed linear functions from [FT22]
pub struct HSNPLinSecretKey<E: PairingEngine> {
    /// the first element of the signing key
    pub a: E::Fr,
    /// the second element of the signing key
    pub b: E::Fr,
    /// the seed of the PRF for the signing key
    pub kappa: E::Fr,
}

///Implements the signature of the HSNP scheme for committed linear functions from [FT22]
pub struct HSNPLinSignature<E: PairingEngine> {
    /// the signed message
    pub msg: E::Fr,
    /// the group element of the signature
    pub lambda: E::G1Affine,
    /// the chameleon hash randomness of the signature
    pub r: E::Fr,
}


///Implements the proof of the HSNP scheme for committed linear functions from [FT22]
pub struct HSNPLinProof<E: PairingEngine> {
    /// the group element of the signature
    pub lambda: E::G1Affine,
    /// the chameleon hash randomness of the signature
    pub r: E::Fr,
    ///the proof of knowledge of the commitment to the result
    pub pok_z: hsnp_pc::PoKProof<E>,
}

impl<E: PairingEngine> ToBytes for HSNPLinProof<E> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> ark_std::io::Result<()> {
        self.lambda.write(&mut writer)?;
        self.r.write(&mut writer)?;
        self.pok_z.write(&mut writer)
    }
}