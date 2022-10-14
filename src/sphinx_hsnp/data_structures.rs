use ark_ec::PairingEngine;
use ark_poly::{GeneralEvaluationDomain, univariate::DensePolynomial as DensePoly};
use ark_poly_commit::hsnp_pc;
use crate::lhs_ped;
use ark_ff::ToBytes;
use ark_std::io::Write;


///Implements the public key of the HSNP scheme [FT22]
pub struct R1csHsnpPublicKey<E: PairingEngine> {
    ///public parameters of the polynomial commitment
    pub pp: hsnp_pc::UniversalParams<E>,
    ///public key of the HSNPLin scheme
    pub pklin: lhs_ped::HSNPLinPublicKey<E>,
}

///Implements the secret key of the HSNP scheme from [FT22]
pub struct R1csHsnpSecretKey<E: PairingEngine> {
    ///secret key of the HSNPLin scheme
    pub sklin: lhs_ped::HSNPLinSecretKey<E>,
}

///Implements the preprocessed public key of the HSNP scheme from [FT22]
pub struct R1csHsnpEvKeyR<'a, E: PairingEngine> {
    /*
    ///public parameters of the polynomial commitment
    pub pp: hsnp_pc::UniversalParams<E>,
    */
    ///hsnp_pc commitment key
    pub ck: hsnp_pc::Powers<'a, E>,
    ///Lagrange-basis commitment key
    pub ck_t: Vec<E::G1Affine>,
    ///z_t polynomial
    pub z_t: DensePoly<E::Fr>,
    ///hsnp_pc verification key
    pub cvk: hsnp_pc::VerifierKey<E>,
    ///public key of the HSNPLin scheme
    pub pklin: lhs_ped::HSNPLinPublicKey<E>,
    ///evaluation domain of size num_constraints
    pub domain_h: GeneralEvaluationDomain<E::Fr>,
    ///preocomputed evaluations of the z_T polynomial in a coset to be used in the division in CPsvec
    pub z_t_evals: Vec<E::Fr>,
}

///Implements the preprocessed verification key of the HSNP scheme from [FT22]
pub struct R1csHsnpVerKeyR<E: PairingEngine> {
    ///public key of the HSNPLin scheme
    pub pklin: lhs_ped::HSNPLinPublicKey<E>,
    ///specialized verification key
    pub vk_t: hsnp_pc::SpecializedVerifierKey<E>,
    ///evaluation domain of size num_constraints
    pub domain_h: GeneralEvaluationDomain<E::Fr>,
}

/// `R1csHsnpSignature` is the signature the HSNP scheme.
pub type R1csHsnpSignature<E> = lhs_ped::HSNPLinSignature<E>;

///Implements the proof of the HSNP scheme from [FT22]
pub struct R1csHsnpProof<E: PairingEngine> {
    ///commitment to the Marlin's witness
    /// (this is for the moment a simulation of Marlin's proof)
    /// to be replaced by Marlin's proof when integrated
    pub com_wprime: hsnp_pc::Commitment<E>,
    ///commitment to the input vector x
    pub com_x: hsnp_pc::Commitment<E>,
    ///CPsvec proof
    pub pi_svec: CPsvecProof<E>,
    ///commitment to the sum of x_i
    pub com_sum: hsnp_pc::Commitment<E>,
    ///HSLin proof
    pub sigproof_sum: lhs_ped::HSNPLinProof<E>,
    ///commitment to z
    pub com_z: hsnp_pc::Commitment<E>,
    ///committed evaluation proof that z=<x, Lag(r)
    pub pi_ev: hsnp_pc::ComEvProof<E>,
    ///HSLin proof that z=<x, Lag(r)>
    pub sigproof_z: lhs_ped::HSNPLinProof<E>,
}

///Implements the proof of the CPsvec CPSNARK from [FT22]
pub struct CPsvecProof<E: PairingEngine> {
    ///commitment to the shifted version of w
    pub com_w: hsnp_pc::Commitment<E>,
    ///commitment to the randomness x
    pub com_0: hsnp_pc::Commitment<E>,
    ///opening proof
    pub pi_0: hsnp_pc::Proof<E>,
    
}

impl<E: PairingEngine> ToBytes for CPsvecProof<E> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> ark_std::io::Result<()> {
        self.com_w.write(&mut writer)?;
        self.com_0.write(&mut writer)?;
        self.pi_0.write(&mut writer)
    }
}