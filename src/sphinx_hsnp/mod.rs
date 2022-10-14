use crate::{Error,lhs_ped};
use ark_std::marker::PhantomData;
use ark_ec::{PairingEngine, AffineCurve, ProjectiveCurve};
use ark_ff::{FftField, FromBytes, to_bytes};
use ark_poly::{univariate::DensePolynomial, UVPolynomial, GeneralEvaluationDomain, EvaluationDomain, Evaluations as EvaluationsOnDomain};
use ark_std::ops::{Sub, Mul, Div};
use ark_poly_commit::hsnp_pc;
use ark_poly_commit::PCUniversalParams;
use ark_std::rand::{RngCore, SeedableRng};
use ark_relations::r1cs::SynthesisError;
use ark_std::*;
use rand_chacha::ChaChaRng;
use blake2::{Blake2s, Digest};

use ark_relations::r1cs::ConstraintSynthesizer;
use ark_marlin::ahp::AHPForR1CS;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

mod data_structures;
pub use data_structures::*;

///Implements the HSNP for R1CS obtained by instantiating the generic HSNP construction of [FT22]
pub struct R1csHsnp<E: PairingEngine, P: UVPolynomial<E::Fr>> {
    _engine: PhantomData<E>,
    _poly: PhantomData<P>,
}

impl<E, P> R1csHsnp<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    for<'a, 'b> &'a P: Sub<&'b P, Output = P>,
{
    /// protocol name for the scheme, to be used to seed FiatShamir rng
    pub const HSNP_NAME: &'static [u8] = b"HSNP_R1CS";

    ///protocol name for CPsvec, to be used to seed FiatShamir rng
    pub const CPSVEC_NAME: &'static [u8] = b"CPsvec";
    
    ///Setup algorithm
    pub fn setup<R: RngCore>(
        max_nconstraints: usize,
        max_nvariables: usize,
        max_non_zero: usize,
        max_inputs: usize,
        rng: &mut R,
    ) -> Result<(R1csHsnpPublicKey<E>, R1csHsnpSecretKey<E>), Error> {

        let max_degree = AHPForR1CS::<E::Fr>::max_degree(max_nconstraints, max_nvariables, max_non_zero).unwrap();

        if max_inputs > max_nconstraints {
            return Err(Error::TooFewConstraints);
        }
        
        let pp = hsnp_pc::HsnpPc::setup(max_degree, true, rng).unwrap();
        
        //extracts the elements (g, gamma_g, h) that are sufficient to
        //commit to a scalar and verify proofs 
        let ck = hsnp_pc::HsnpPc::get_scalar_comkey(&pp).unwrap();
        
        //generates the universal relation for Pedersen-committed outputs of
        //linear functions. The relation is simply defined by the commitment
        //key of the Pedersen commitment
        let hslin_univrel = lhs_ped::UniversalPedLinRelation {
            ck,
        };
        
        //generates the public and secret key of the HSNPLin scheme
        let (pklin, sklin) = lhs_ped::HSNPLin::setup(&hslin_univrel,rng);
        
        let pk = R1csHsnpPublicKey {
            pp,
            pklin,
        };

        let sk = R1csHsnpSecretKey {
            sklin,
        };
        
        Ok((pk, sk))
    }

    ///Outputs a preprocessed evaluation key and verification key for a given relation
    pub fn preprocess<C: ConstraintSynthesizer<E::Fr>>(
        pk: &R1csHsnpPublicKey<E>,
        c: C,
        //num_constraints: usize,
        num_inputs: usize,
    ) -> Result<(R1csHsnpEvKeyR<E>, R1csHsnpVerKeyR<E>), Error> {

        let max_hiding_bound = num_inputs+2;

        let index = AHPForR1CS::index(c).unwrap();
        if pk.pp.max_degree() < index.max_degree() {
            return Err(Error::IndexTooLarge);
        }

        let supported_degree = index.max_degree(); //num_constraints;
        
        //we define the number of constraints, which is actually the max R1CS dimension in Marlin
        let num_constraints = cmp::max(index.index_info.num_constraints, index.index_info.num_variables);
        
        
        let (ck, vk) = hsnp_pc::HsnpPc::trim(&pk.pp, supported_degree, max_hiding_bound).unwrap();
        let (z_t, vk_t) = hsnp_pc::HsnpPc::specialize(&pk.pp, &vk, num_constraints, num_inputs).unwrap();
        
        //generate FFT interpolation domain of size |H|=num_constraints
        //(this is the same interpolation domain of Marlin)
        let domain_h = GeneralEvaluationDomain::<E::Fr>::new(num_constraints)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge).unwrap();
        
        //generate Lagrange-basis commitment key
        //deriving preprocessed commitment key to be the first t Lagrange polynomials in the exponent
        let generator = E::Fr::get_root_of_unity(domain_h.size()).unwrap();
        let mut powers_of_g_proj = Vec::new();
        for k in 0..domain_h.size() {
            powers_of_g_proj.push(pk.pp.powers_of_g[k].into_projective());
        }
        let proj_ck_t = hsnp_pc::HsnpPc::<E,P>::compute_lagrange_comkey(&powers_of_g_proj, generator, domain_h.size())[..num_inputs].to_vec();
        let ck_t = E::G1Projective::batch_normalization_into_affine(&proj_ck_t);
        
        //precompute evaluations of the z_T polynomial to be used in the division in CPsvec
        let z_t_evals = domain_h.coset_fft(&z_t);
        
        let ek_r = R1csHsnpEvKeyR{
            //pp: pk.pp.clone(),
            ck,
            ck_t,
            z_t,
            cvk: vk,
            pklin: pk.pklin.clone(),
            domain_h,
            z_t_evals,
        };

        let vk_r = R1csHsnpVerKeyR{
            pklin: pk.pklin.clone(),
            vk_t, //this includes vk
            domain_h,
        };

        Ok((ek_r, vk_r))
    }

    ///Outputs a preprocessed evaluation key and verification key for a specific t
    pub fn specialize_to_t<'a, C: ConstraintSynthesizer<E::Fr>>(
        pk: &'a R1csHsnpPublicKey<E>,
        c: C,
        ck_maxt: &'a [E::G1Affine],
        num_inputs: usize,
    ) -> Result<(R1csHsnpEvKeyR<'a, E>, R1csHsnpVerKeyR<E>), Error> {

        let max_hiding_bound = num_inputs+2;

        let index = AHPForR1CS::index(c).unwrap();
        if pk.pp.max_degree() < index.max_degree() {
            return Err(Error::IndexTooLarge);
        }

        let supported_degree = index.max_degree(); //num_constraints;
        
        //we define the number of constraints, which is actually the max R1CS dimension in Marlin
        let num_constraints = cmp::max(index.index_info.num_constraints, index.index_info.num_variables);
        
        
        let (ck, vk) = hsnp_pc::HsnpPc::trim(&pk.pp, supported_degree, max_hiding_bound).unwrap();
        let (z_t, vk_t) = hsnp_pc::HsnpPc::specialize(&pk.pp, &vk, num_constraints, num_inputs).unwrap();
        
        //generate FFT interpolation domain of size |H|=num_constraints
        //(this is the same interpolation domain of Marlin)
        let domain_h = GeneralEvaluationDomain::<E::Fr>::new(num_constraints)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge).unwrap();
        
        //let proj_ck_t = hsnp_pc::HsnpPc::<E,P>::compute_lagrange_comkey(&powers_of_g_proj, generator, domain_h.size())[..num_inputs].to_vec();
        //let ck_t = E::G1Projective::batch_normalization_into_affine(&proj_ck_t);

        let ck_t = ck_maxt[..num_inputs].to_vec();
        
        
        //precompute evaluations of the z_T polynomial to be used in the division in CPsvec
        let z_t_evals = domain_h.coset_fft(&z_t);
        
        let ek_r = R1csHsnpEvKeyR{
            //pp: pk.pp.clone(),
            ck,
            ck_t,
            z_t,
            cvk: vk,
            pklin: pk.pklin.clone(),
            domain_h,
            z_t_evals,
        };

        let vk_r = R1csHsnpVerKeyR{
            pklin: pk.pklin.clone(),
            vk_t, //this includes vk
            domain_h,
        };

        Ok((ek_r, vk_r))
    }

    

    ///Sign algorithm: it simply invokes the sign of the HSNPLin scheme
    pub fn sign(
        pk: &R1csHsnpPublicKey<E>,
        sk: &R1csHsnpSecretKey<E>,
        lab: u128,
        msg: E::Fr,
    ) -> R1csHsnpSignature<E> {
        let sig = lhs_ped::HSNPLin::sign(&pk.pklin, &sk.sklin, lab, msg);
        sig
    }

    ///Evaluation algorithm
    /// For now it does not take an actual relation but only simulates all the steps of the proof but the Marlin's one
    pub fn eval<R: RngCore>(
        ek_r: &R1csHsnpEvKeyR<E>,
        labs: &[u128],
        y: &[E::Fr],
        sigs: &[R1csHsnpSignature<E>],
        w: &[E::Fr],
        rng: &mut R,
    ) -> R1csHsnpProof<E> {

        //============= (Step 1) =============//
        //Extract the vector x of signed inputs and commit to x

        let x = sigs.iter().map(|x| x.msg).collect::<Vec<_>>();
        let x_pol = EvaluationsOnDomain::from_vec_and_domain(x.clone(), ek_r.domain_h).interpolate();
        let (com_x, opn_x) = hsnp_pc::HsnpPc::<E,DensePolynomial<E::Fr>>::commit_from_lagrange_representation(&ek_r.ck, &ek_r.ck_t, &x, Some(0), Some(rng)).unwrap();
        
        //============= (END Step 1: Com(x)) =============//


        //============= (Step 2) =============//
        //Generate CP_R proof of  (y, x, w) \in R
        //Here we explicitly write our CP_R based on Marlin+CPsvec
        
            //===Generate Marlin's proof of (y, x, w) \in R
            //For now this is only simulated by generating the commitment to w'=(x,w)

            //generate w'=(x,w)
        
            
            let w_prime_vec = [x.clone(), w.to_vec()].concat();

            //interpolate w' -> w'(X)
            let w_prime = EvaluationsOnDomain::from_vec_and_domain(w_prime_vec, ek_r.domain_h).interpolate();
            
            //commit to w'
            let (com_wprime, opn_wprime) = hsnp_pc::HsnpPc::<E,DensePolynomial<E::Fr>>::commit(&ek_r.ck, &w_prime, Some(2), Some(rng)).unwrap();
            
            //===END Marlin's proof

            //generate CPsvec proof
            let div_time = start_timer!(|| "Computing division");
            //evaluate x and wprime on domain's coset of size H
            let x_evals = ek_r.domain_h.coset_fft(&x_pol);
            let mut w_evals = ek_r.domain_h.coset_fft(&w_prime);
            w_evals = cfg_into_iter!(0..ek_r.domain_h.size())
                .map(|k| {
                    (w_evals[k] - &x_evals[k])/&ek_r.z_t_evals[k]
                }).collect();
            //compute w by interpolation in the coset
            ek_r.domain_h.coset_ifft_in_place(&mut w_evals);
            let w_pol = DensePolynomial::<E::Fr>::from_coefficients_vec(w_evals);
            end_timer!(div_time);
            
        
            //-- Com(w) hiding bound 2
            let (com_w, opn_w) = hsnp_pc::HsnpPc::<E,DensePolynomial<E::Fr>>::commit(&ek_r.ck, &w_pol, Some(2), Some(rng)).unwrap();

            // com_0 = h1^r_0 s.t. r_0(X) = r_w'(X) - r_x - r_w(X)*Z_t(X)
            //compute r_0
            let rand_0 = &(&opn_wprime.blinding_polynomial - &opn_x.blinding_polynomial) - &(opn_w.blinding_polynomial.mul(&ek_r.z_t));
            //compute h_1^r_0
            let (com_0, opn_0) = hsnp_pc::HsnpPc::<E,DensePolynomial<E::Fr>>::deterministic_commit(&ek_r.ck, &DensePolynomial::<E::Fr>::zero(), &rand_0).unwrap();
        
            let bytes = to_bytes![&Self::CPSVEC_NAME, &com_x.0, &com_w.0, &com_0.0, &com_wprime.0].unwrap();
            let rho = Self::fiat_shamir_hashfn(&bytes);
            
            let pi_0 = hsnp_pc::HsnpPc::<E,DensePolynomial<E::Fr>>::open(&ek_r.ck, &DensePolynomial::<E::Fr>::zero(), rho, &opn_0).unwrap();

            let pi_svec = CPsvecProof {
                com_w,
                com_0,
                pi_0,
            };

            //===END CPsvec proof
        //============= (END Step 2: CP_R proof) =============//


        //============= (Step 3) =============//
        // Generate HSNPLin proof for the sum_i x_i
        // This is a way to prove possession of signatures for each labeled input

        let ones = vec![E::Fr::one(); sigs.len()];
        let x_sum: E::Fr = x.iter().sum();

        let linrel_sum = lhs_ped::PedLinRelation::<E> {
            signed_input_size: sigs.len(),
            f: ones,
        };
        
        let x_sum_pol = DensePolynomial::<E::Fr>::from_coefficients_vec(vec![x_sum]);
        let (com_sum, opn_sum) = hsnp_pc::HsnpPc::<E,DensePolynomial<E::Fr>>::commit(&ek_r.ck, &x_sum_pol, Some(0), Some(rng)).unwrap();
        let sigproof_sum = lhs_ped::HSNPLin::<E, DensePolynomial<E::Fr>>::eval(&ek_r.pklin, &linrel_sum, &com_sum, &x_sum_pol, &opn_sum, &sigs, rng);
        
        //============= (END Step 3: HSNPLin sum proof) =============//


        //============= (Step 4) =============//
        // Generate random oracle challenge
        
        //TODO: after integration with Marlin add pi_y to the Fiat-Shamir challenge

        let bytes = to_bytes![&Self::HSNP_NAME, y, &com_x.0, &pi_svec, labs, &com_sum.0, &sigproof_sum].unwrap();
        let r = Self::fiat_shamir_hashfn(&bytes);
        
        //============= (END Step 4: r) =============//


        //============= (Step 5) =============//
        // Compute z = <x, Lag(r)> and its commitment com_z = Com(z)
        
        let s = eval_lagrange(ek_r.domain_h, x.len(), r);
        let z: E::Fr = s.iter().zip(x.iter()).map(|(a,b)| *a*b).sum();
        let z_pol = DensePolynomial::<E::Fr>::from_coefficients_vec(vec![z]);
        let (com_z, opn_z) = hsnp_pc::HsnpPc::<E,DensePolynomial<E::Fr>>::commit(&ek_r.ck, &z_pol, Some(0), Some(rng)).unwrap();
        
        //============= (END Step 5: Com(z)) =============//


        //============= (Step 6) =============//
        // Generate CPev proof
        
        let pi_ev = hsnp_pc::HsnpPc::<E,DensePolynomial<E::Fr>>::committed_evaluation_prove(&ek_r.ck, &ek_r.cvk, r, &com_x, &x_pol, &opn_x, &com_z, &z_pol, &opn_z, rng).unwrap();

        //============= (END Step 6: CPev) =============//
        

        //============= (Step 7) =============//
        // Generate HSNPLin proof for com_z = Com(<x,Lag(r)>)
        
        let linrel_s = lhs_ped::PedLinRelation::<E> {
            signed_input_size: sigs.len(),
            f: s,
        };
        let sigproof_z = lhs_ped::HSNPLin::<E, DensePolynomial<E::Fr>>::eval(&ek_r.pklin, &linrel_s, &com_z, &z_pol, &opn_z, &sigs, rng);
        
        //============= (END Step 7: HSNPLin z) =============//
        

        //Set proof
        let proof = R1csHsnpProof {
            com_wprime,
            com_x,
            pi_svec,
            com_sum,
            sigproof_sum,
            com_z,
            pi_ev,
            sigproof_z,
        };
        proof
    }

    /// Verification algorithm
    pub fn verify(
        vk_r: &R1csHsnpVerKeyR<E>,
        labs: &Vec<u128>,
        y: &[E::Fr],
        sigproof: &R1csHsnpProof<E>,
    ) -> bool {
    
        // == check CPsvec proof ==
        
        //   e(com_w' / (com_x com_0), h) = e(com_w, z_t)
        let inner = sigproof.com_wprime.0.into_projective() - &sigproof.com_x.0.into_projective() - &sigproof.pi_svec.com_0.0.into_projective();
        let lhs = E::pairing(inner, vk_r.vk_t.vk.h);
        let rhs = E::pairing(sigproof.pi_svec.com_w.0, vk_r.vk_t.com_zt);
        let svec_ver1 = lhs == rhs;
        
        //check that com_0 commits to 0 by checking an evaluation in random point rho
        let bytes = to_bytes![&Self::CPSVEC_NAME, &sigproof.com_x.0, &sigproof.pi_svec.com_w.0, &sigproof.pi_svec.com_0.0, &sigproof.com_wprime.0].unwrap();
        let rho = Self::fiat_shamir_hashfn(&bytes);
        let svec_ver2 = hsnp_pc::HsnpPc::check(&vk_r.vk_t.vk, &sigproof.pi_svec.com_0, rho, E::Fr::zero(), &sigproof.pi_svec.pi_0).unwrap();
        
        let ver_svec = svec_ver1 && svec_ver2;

        // == check HSLin proof for sum ==
        
        //defines the sum linear function as vector of 1s
        let ones = vec![E::Fr::one(); labs.len()];
        let linrel_sum = lhs_ped::PedLinRelation::<E> {
            signed_input_size: labs.len(),
            f: ones,
        };
        //verify that com_sum commits to the sum of the signed values
        // (recall: this is a way to check that all the labels have been signed)
        let ver_hslin_sum = lhs_ped::HSNPLin::ver(&vk_r.pklin, &linrel_sum, labs, &sigproof.com_sum, &sigproof.sigproof_sum);
        
        // == compute random oracle challenge ==
        let bytes = to_bytes![&Self::HSNP_NAME, y, &sigproof.com_x.0, &sigproof.pi_svec, labs, &sigproof.com_sum.0, &sigproof.sigproof_sum].unwrap();
        let r = Self::fiat_shamir_hashfn(&bytes);
        
        // == check CPev proof ==
        let ver_ev = hsnp_pc::HsnpPc::committed_evaluation_ver(&vk_r.vk_t.vk, r, &sigproof.com_x, &sigproof.com_z, &sigproof.pi_ev).unwrap();
        
        //check HSLin proof for z
        
        let s = eval_lagrange(vk_r.domain_h, labs.len(), r);
        let linrel_s = lhs_ped::PedLinRelation::<E> {
            signed_input_size: labs.len(),
            f: s,
        };
        let ver_hslin_z = lhs_ped::HSNPLin::ver(&vk_r.pklin, &linrel_s, labs, &sigproof.com_z, &sigproof.sigproof_z);
    
        //accept iff all the above verifications succeeded
        ver_svec && ver_hslin_sum && ver_ev && ver_hslin_z
    }

    ///precompute the vector ( R_lab=H(lab) )
    pub fn precompute_r_labs(
        labs: &Vec<u128>,
    ) -> Vec<E::G1Affine> {
        let r_labs = lhs_ped::HSNPLin::<E, P>::precompute_r_labs(labs);
        r_labs
    }

    /// Verification algorithm
    pub fn verify_with_precomputed_r_labs(
        vk_r: &R1csHsnpVerKeyR<E>,
        labs: &Vec<u128>,
        r_labs: &Vec<E::G1Affine>,
        y: &[E::Fr],
        sigproof: &R1csHsnpProof<E>,
    ) -> bool {
    
        // == check CPsvec proof ==
        
        //   e(com_w' / (com_x com_0), h) = e(com_w, z_t)
        let inner = sigproof.com_wprime.0.into_projective() - &sigproof.com_x.0.into_projective() - &sigproof.pi_svec.com_0.0.into_projective();
        let lhs = E::pairing(inner, vk_r.vk_t.vk.h);
        let rhs = E::pairing(sigproof.pi_svec.com_w.0, vk_r.vk_t.com_zt);
        let svec_ver1 = lhs == rhs;
        
        //check that com_0 commits to 0 by checking an evaluation in random point rho
        let bytes = to_bytes![&Self::CPSVEC_NAME, &sigproof.com_x.0, &sigproof.pi_svec.com_w.0, &sigproof.pi_svec.com_0.0, &sigproof.com_wprime.0].unwrap();
        let rho = Self::fiat_shamir_hashfn(&bytes);
        let svec_ver2 = hsnp_pc::HsnpPc::check(&vk_r.vk_t.vk, &sigproof.pi_svec.com_0, rho, E::Fr::zero(), &sigproof.pi_svec.pi_0).unwrap();
        
        let ver_svec = svec_ver1 && svec_ver2;

        // == check HSLin proof for sum ==
        
        //defines the sum linear function as vector of 1s
        let ones = vec![E::Fr::one(); labs.len()];
        let linrel_sum = lhs_ped::PedLinRelation::<E> {
            signed_input_size: labs.len(),
            f: ones,
        };
        //verify that com_sum commits to the sum of the signed values
        // (recall: this is a way to check that all the labels have been signed)
        let ver_hslin_sum = lhs_ped::HSNPLin::ver_with_precomputed_r_labs(&vk_r.pklin, &linrel_sum, r_labs, &sigproof.com_sum, &sigproof.sigproof_sum);
        
        // == compute random oracle challenge ==
        let bytes = to_bytes![&Self::HSNP_NAME, y, &sigproof.com_x.0, &sigproof.pi_svec, labs, &sigproof.com_sum.0, &sigproof.sigproof_sum].unwrap();
        let r = Self::fiat_shamir_hashfn(&bytes);
        
        // == check CPev proof ==
        let ver_ev = hsnp_pc::HsnpPc::committed_evaluation_ver(&vk_r.vk_t.vk, r, &sigproof.com_x, &sigproof.com_z, &sigproof.pi_ev).unwrap();
        
        //check HSLin proof for z
        
        let s = eval_lagrange(vk_r.domain_h, labs.len(), r);
        let linrel_s = lhs_ped::PedLinRelation::<E> {
            signed_input_size: labs.len(),
            f: s,
        };
        let ver_hslin_z = lhs_ped::HSNPLin::ver_with_precomputed_r_labs(&vk_r.pklin, &linrel_s, r_labs, &sigproof.com_z, &sigproof.sigproof_z);
    
        //accept iff all the above verifications succeeded
        ver_svec && ver_hslin_sum && ver_ev && ver_hslin_z
    }

    ///hash function to generate random oracle challenge in Fr
    pub(crate) fn fiat_shamir_hashfn(
        bytes: &[u8],
    ) -> E::Fr {
        let chg = Blake2s::digest(bytes);
    
        let seed: [u8; 32] = FromBytes::read(&*chg).expect("failed to get [u32; 8]");
        let mut hrng = ChaChaRng::from_seed(seed);
        let rho = E::Fr::rand(&mut hrng);
        rho
    }
}

/// Evaluates a vector of the first t Lagrange polynomials 
/// of degree n-1 in a point r,
/// with running time O(max(t, log n))
fn eval_lagrange<F: FftField>(domain: GeneralEvaluationDomain<F>, t: usize, r: F) -> Vec<F> {

    use ark_ff::fields::batch_inversion;
    
    let size = domain.size() as u64;
    //domain size as field element
    let size_as_field_elm = F::from(size);
    
    //compute vanishing polynomial at r: v_H(r)
    let van_at_r = domain.evaluate_vanishing_polynomial(r);
    
    //set generator and first element of the domain
    let gen = F::get_root_of_unity(domain.size()).unwrap();
    let gen_inv = gen.pow(&[size -1]);
    let domain_offset = F::one();
    
    //coefficient^-1 = (m(r - g^i)) * g^-i v_H(r)^-1

    let mut s = vec![F::zero(); t];
    let mut h_i = domain_offset;
    let mut v_0_inv = size_as_field_elm * van_at_r.inverse().unwrap();
    for i in 0..t {
        //set h_i = g^i
        if i>0 {
            h_i = h_i * gen;
            v_0_inv *= gen_inv;
        }
        //denominator = m (r - h_i)
        let denom_i =  r - h_i;
        //coefficient = h_i v_H(r) / (m(r - h_i))
        s[i] = v_0_inv * denom_i;
        //s.push(coeff_i);
    }
    batch_inversion(s.as_mut_slice());


    s

}