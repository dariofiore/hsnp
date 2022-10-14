
use ark_ec::{PairingEngine, AffineCurve, ProjectiveCurve, msm::VariableBaseMSM};
use ark_ff::{One, UniformRand, PrimeField};
use ark_std::marker::PhantomData;
use ark_std::rand::{RngCore, SeedableRng};
use rand_chacha::ChaChaRng;
use ark_ff::to_bytes;
use ark_ff::FromBytes;
use blake2::{Blake2s, Digest};
use ark_poly::UVPolynomial;
use ark_poly_commit::hsnp_pc;
use ark_std::ops::Div;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

mod data_structures;
pub use data_structures::*;

///Implements the HSNP for committed linear functions from [FT22]
pub struct HSNPLin<E: PairingEngine, P: UVPolynomial<E::Fr>> {
    _engine: PhantomData<E>,
    _poly: PhantomData<P>,
}

impl<E, P> HSNPLin<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    ///protocol name to be used to seed the hash function 
    /// that publicly generates group elements for signing
    /// and verifying
    pub const HSNPLIN_HASH: &'static [u8] = b"HSNPLin-hashseed";

    /// name to be used to seed the PRF that generates the signatures
    /// randomness from a secret seed
    pub const HSNPLIN_PRF: &'static [u8] = b"HSNPLin-prfseed";

        
    ///Setup algorithm of the scheme
    pub fn setup<R: RngCore>(
        univ_rel: &UniversalPedLinRelation<E>,
        rng: &mut R,
    ) -> (HSNPLinPublicKey<E>, HSNPLinSecretKey<E>) {

        let (a, b) = (E::Fr::rand(rng), E::Fr::rand(rng));
        let kappa = E::Fr::rand(rng);
        
        let b_g = univ_rel.ck.g.mul(b).into_affine();
        let gamma1 = univ_rel.ck.gamma_g.mul(a).into_affine();
        let a_inv = E::Fr::one() / &a;
        let gamma2 = univ_rel.ck.h.mul(a_inv).into_affine();

        let pk = HSNPLinPublicKey {
            ck: univ_rel.ck.clone(),
            gamma1,
            gamma2,
            b_g,
        };

        let sk = HSNPLinSecretKey {
            a,
            b,
            kappa,
        };
        (pk, sk)
    }

    ///Signing algorithm of the scheme
    pub fn sign(
        pk: &HSNPLinPublicKey<E>,
        sk: &HSNPLinSecretKey<E>,
        lab: u128,
        msg: E::Fr,
    ) -> HSNPLinSignature<E> {
        
        //R_lab = H(lab)
        let r_lab = Self::lab_hash(lab);

        //r = PRF_k(lab)
        let r = Self::prf(sk.kappa, lab);

        //compute lambda = (R_lab g_1^{m + rb})^a
        let mut exponent = msg + sk.b * r;
        exponent *= sk.a;
        let lambda = pk.ck.g.mul(exponent) + &r_lab.mul(sk.a);
        
        let sig = HSNPLinSignature {
            msg,
            lambda: lambda.into_affine(),
            r,
        };
        sig
    }
    
    ///Evaluation algorithm
    pub fn eval<R: RngCore>(
        pk: &HSNPLinPublicKey<E>,
        rel: &PedLinRelation<E>,
        com_z: &hsnp_pc::Commitment<E>,
        z: &P,
        rand_z: &hsnp_pc::Randomness<E::Fr, P>,
        sigs: &[HSNPLinSignature<E>],
        rng: &mut R,
    ) -> HSNPLinProof<E> {

        //compute a ZKPoK
        let ck = hsnp_pc::ScalarComKey {
            g: pk.ck.g,
            gamma_g: pk.ck.gamma_g,
            h: pk.ck.h,
        };
        let pi_z = hsnp_pc::HsnpPc::pok_scalar_prove(&ck, com_z, z, rand_z, rng).unwrap();
        
        let opn_z = rand_z.blinding_polynomial.coeffs()[0];

        //performs multi-exponentiation for Lambda

        let int_coeffs = rel.f.iter()
            .map(|s| s.into_repr())
                .collect::<Vec<_>>();
        
        let mut lambda = VariableBaseMSM::multi_scalar_mul(
            &sigs.iter().map(|x| x.lambda).collect::<Vec<_>>(),
            &int_coeffs,
        );

        let lambda_z = pk.gamma1.mul(opn_z);
        lambda = lambda + &lambda_z;
        //computes <r,s>
        let r: E::Fr = sigs.iter().zip(rel.f.iter()).map(|(x,y)| x.r*y).sum();
        
        let proof = HSNPLinProof {
            lambda: lambda.into_affine(),
            r,
            pok_z: pi_z,
        };
        proof
    }

    ///Verification algorithm
    pub fn ver(
        pk: &HSNPLinPublicKey<E>,
        rel: &PedLinRelation<E>,
        labs: &Vec<u128>,
        com_z: &hsnp_pc::Commitment<E>,
        proof: &HSNPLinProof<E>,
    ) -> bool {

        //compute prod_lab R_lab^s_lab with R_lab = H(lab)
        //this is slow due to the hashing based on rand() generation
        let r_labs_comp_time = start_timer!(|| "Computing {H(lab)}");
        //the following initialization is to fix the vector's length
        let mut r_labs = Vec::<E::G1Affine>::with_capacity(labs.len());
        //the following is a parallel-friendly iterator
        r_labs = cfg_into_iter!(0..labs.len()).map(|k| { Self::lab_hash(labs[k]) }).collect();
        end_timer!(r_labs_comp_time);
        let r_labs_prod_time = start_timer!(|| "Computing prod_i H(lab_i)^f_i");
        let r_f = VariableBaseMSM::multi_scalar_mul(
            &r_labs,
            &rel.f.iter().map(|s| s.into_repr()).collect::<Vec<_>>(),
        ).into_affine();
        end_timer!(r_labs_prod_time);
        
        let lhs = E::pairing(proof.lambda, pk.gamma2);
        let inner = r_f.into_projective() + &com_z.0.into_projective() + &pk.b_g.mul(proof.r);
        let rhs = E::pairing(inner.into_affine(), pk.ck.h);

        let ck = hsnp_pc::ScalarComKey {
            g: pk.ck.g,
            gamma_g: pk.ck.gamma_g,
            h: pk.ck.h,
        };
        
        let pok_ver = hsnp_pc::HsnpPc::pok_scalar_ver(&ck, com_z, &proof.pok_z).unwrap();

        lhs == rhs && pok_ver
    }

    ///Verification algorithm
    pub fn ver_with_precomputed_r_labs(
        pk: &HSNPLinPublicKey<E>,
        rel: &PedLinRelation<E>,
        r_labs: &Vec<E::G1Affine>,
        com_z: &hsnp_pc::Commitment<E>,
        proof: &HSNPLinProof<E>,
    ) -> bool {

        let r_labs_prod_time = start_timer!(|| "Computing prod_i H(lab_i)^f_i");
        let r_f = VariableBaseMSM::multi_scalar_mul(
            &r_labs,
            &rel.f.iter().map(|s| s.into_repr()).collect::<Vec<_>>(),
        ).into_affine();
        end_timer!(r_labs_prod_time);
        
        let lhs = E::pairing(proof.lambda, pk.gamma2);
        let inner = r_f.into_projective() + &com_z.0.into_projective() + &pk.b_g.mul(proof.r);
        let rhs = E::pairing(inner.into_affine(), pk.ck.h);

        let ck = hsnp_pc::ScalarComKey {
            g: pk.ck.g,
            gamma_g: pk.ck.gamma_g,
            h: pk.ck.h,
        };
        
        let pok_ver = hsnp_pc::HsnpPc::pok_scalar_ver(&ck, com_z, &proof.pok_z).unwrap();

        lhs == rhs && pok_ver
    }
    

    ///this is a hash to G1
    pub(crate) fn lab_hash(
        lab: u128,
    ) -> E::G1Affine {
        let name = Self::HSNPLIN_HASH;
        let mut ctr: u64 = 0;
        let mut g_elm = None;
        while g_elm == None {
            ctr += 1;
            let bytes = to_bytes![name, lab, ctr].unwrap();
            let bytehash = Blake2s::digest(&bytes);
            g_elm = E::G1Affine::from_random_bytes(&bytehash);
            
        }
        g_elm.unwrap()

    }

    ///precompute the vector ( R_lab=H(lab) )
    pub fn precompute_r_labs(
        labs: &Vec<u128>,
    ) -> Vec<E::G1Affine> {
        let mut r_labs = Vec::<E::G1Affine>::with_capacity(labs.len());
        //the following is a parallel-friendly iterator
        r_labs = cfg_into_iter!(0..labs.len()).map(|k| { Self::lab_hash(labs[k]) }).collect();
        r_labs
    }

    ///this is a PRF that outputs a pseudorandom field element
    pub(crate) fn prf(
        seed: E::Fr,
        lab: u128,
    ) -> E::Fr {
        let name = Self::HSNPLIN_PRF;
        let bytes = to_bytes![name, seed, lab].unwrap();
        let bytehash = Blake2s::digest(&bytes);
        let seed: [u8; 32] = FromBytes::read(&*bytehash).expect("failed to get [u32; 8]");
        let mut hrng = ChaChaRng::from_seed(seed);
        let r = E::Fr::rand(&mut hrng);
        r
    }
}
