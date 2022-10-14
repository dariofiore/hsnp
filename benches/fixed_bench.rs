use ark_bls12_381::*;
use ark_ff::{UniformRand, Zero, PrimeField};
use ark_poly::{univariate::DensePolynomial as DensePoly};
//use ark_poly_commit::hsnp_pc;
//use ark_poly_commit::kzg10::*;
//use ark_ec::msm::{VariableBaseMSM};
//use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_std::*;
//use rayon::prelude::*;
//use ark_relations::r1cs::SynthesisError;
//use ark_ec::PairingEngine;
//use ark_std::ops::{Sub, Mul};
//use blake2::{Blake2s, Digest};

//use hsnp::lhs_ped;
use hsnp::sphinx_hsnp;
use std::env;
use std::str::FromStr;

use ark_relations::{
    lc,
    r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError},
};
use indicatif::{ProgressBar, ProgressStyle};

const NUM_REPETITIONS: usize = 10;

fn main() {
    let args: Vec<String> = env::args().collect();
    /*
    eprintln!(
        "Running HSNP benchmarks with {} repetitions",
        NUM_REPETITIONS
    );*/
    bench_hsnp(args);
}

macro_rules! hsnp_benchmark {
    ($bench_name:ident, $bench_field:ty, $bench_pairing_engine:ty, $bench_precompute_hash:ident, $bench_tests_set:ident) => {
        let rng = &mut ark_std::test_rng();
        type Poly = DensePoly<$bench_field>;
        type HS = sphinx_hsnp::R1csHsnp<$bench_pairing_engine, Poly>;

        let precompute_hash = $bench_precompute_hash;
        let tests_set = $bench_tests_set;
        //compute the max parameters in the tests set
        let mut max_nconstraints = 0;
        let mut max_inputs = 0;
        let mut max_nvariables= 0;
        for test in tests_set.iter() {
            max_inputs = cmp::max(max_inputs, test[0]);
            max_nconstraints = cmp::max(max_nconstraints, test[1]);
            max_nvariables = cmp::max(max_nvariables, test[2]);
        }
        //fix upper bounds for the universal parameters
        let max_non_zero = 3*cmp::max(max_nconstraints, max_nvariables);
        
        eprintln!("[+] Max inputs: {}, max constraints: {}, max variables: {}", max_inputs, max_nconstraints, max_nvariables);

        //generate the public and secret key
        eprint!("[+] Generating Keys...");
        
        let start_srs_gen = ark_std::time::Instant::now();
        let (pk, sk) = HS::setup(max_nconstraints, max_nvariables, max_non_zero, max_inputs, rng).unwrap();
        eprintln!("...done in {} ms.",start_srs_gen.elapsed().as_millis() as u128);
        
        //generate a dummy circuit with `num_constraints` constraints and `num_constraints`-1 variables
        let c = DummyCircuit::<$bench_field> {
            a: Some(<$bench_field>::rand(rng)),
            b: Some(<$bench_field>::rand(rng)),
            num_variables: max_nvariables,
            num_constraints: max_nconstraints,
        };

        eprint!("\t[+] Specializing SRS to t={} inputs and R1CS with dimension ({}x{})...", max_inputs, c.num_constraints, c.num_variables);
        //specialize the keys to the relation for the above circuit
        //and for signed inputs of length `t`
        let start_srs_spec = ark_std::time::Instant::now();
        let (ek_r_maxt, _) = HS::preprocess(&pk, c, max_inputs).unwrap();
        eprintln!("...done in {} ms.",start_srs_spec.elapsed().as_millis() as u128);

        for test in tests_set.iter() {
            eprintln!("[+] Running t={}, n={} ...", test[0], test[1]);

            //fix parameters
            let t = test[0];
            let num_constraints  = test[1];
            let num_variables  = test[2];

            //respecialize
            let (ek_r, vk_r) = HS::specialize_to_t(&pk, c, &ek_r_maxt.ck_t, t).unwrap();

            //initialize time counters
            let mut signing_time = 0 as u128;
            
            //let mut proving_time = 0 as u128;
            //let mut verif_time = 0 as u128;
            let mut proving_times = vec![0 as u128; NUM_REPETITIONS];
            let mut verif_times = vec![0 as u128; NUM_REPETITIONS];
            
            let pb = ProgressBar::new(NUM_REPETITIONS as u64);
            pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:.cyan/blue}] {pos}/{len}")
            .progress_chars("#>-"));
            
            pb.println(format!("\t[+] Signing, proving and verifying iterations"));
            
            //loop for every repetition
            for k in 0..NUM_REPETITIONS {
                
                //sample w and x at random
                let mut x = vec![Fr::zero();t];
                let mut w = vec![Fr::zero(); num_constraints - t];
                x = x.iter().map(|_| {Fr::rand(rng)}).collect();
                w = w.iter().map(|_| {Fr::rand(rng)}).collect();
                
                //generate signatures and labels
                let mut sigs = Vec::with_capacity(t);
                let mut labs = Vec::with_capacity(t);
                let start_signing = ark_std::time::Instant::now();
                for i in 1..=t {
                    let sig = HS::sign(&pk, &sk, (i as u128), x[i-1]);
                    sigs.push(sig);
                    labs.push((i as u128));
                }
                signing_time += start_signing.elapsed().as_micros() as u128;


                //proving
                let y = vec![Fr::rand(rng); 1];
                
                let start_proving = ark_std::time::Instant::now();
                let proof = HS::eval(&ek_r, &labs, &y, &sigs, &w, rng);
                
                proving_times[k] = start_proving.elapsed().as_micros() as u128;
                //proving_time += proving_times[k];

                //verifying
                let r_labs = HS::precompute_r_labs(&labs);
                let start_verifying = ark_std::time::Instant::now();
                let ver = if(precompute_hash) { HS::verify_with_precomputed_r_labs(&vk_r, &labs, &r_labs, &y, &proof) } else { HS::verify(&vk_r, &labs, &y, &proof) };
                //let ver = HS::verify(&vk_r, &labs, &y, &proof);
                //let ver = HS::verify_with_precomputed_r_labs(&vk_r, &labs, &r_labs, &y, &proof);
                verif_times[k] = start_verifying.elapsed().as_micros() as u128;
                //verif_time += verif_times[k];
                if ver == false {
                    eprintln!("verification = {}", ver);
                }
                
                pb.inc(1);
            }
            
            signing_time = signing_time / (NUM_REPETITIONS * t) as u128; //signing time per signature
            
            //proving_time = proving_time / NUM_REPETITIONS as u128;
            //verif_time = verif_time / NUM_REPETITIONS as u128;
            
            //compute the medians
            proving_times.sort();
            let mid = proving_times.len()/2;
            let median_proving_time  = proving_times[mid];

            verif_times.sort();
            let mid = verif_times.len()/2;
            let median_verif_time  = verif_times[mid];
            
            pb.finish_with_message("\t...done");

            println!("{}, {}, {}, {}, {}, {}",
                t,
                num_constraints,
                num_variables,
                signing_time,
                median_proving_time,
                median_verif_time
            );
        }
    };
}

fn bench_hsnp(args: Vec<String>) {
    
    let precompute_hash: bool = if <bool as FromStr>::from_str(&args[1]).is_err(){false} else {FromStr::from_str(&args[1]).unwrap()};
    let benchmark = &args[2];
    let first_t: usize = args[3].parse().unwrap();
    let last_t: usize = args[4].parse().unwrap();
    
    let mut test_set  = Vec::with_capacity(last_t - first_t + 1);
    if benchmark == "var" {
        for j in first_t..=last_t {
            test_set.push([1<<j, (1 << j) + 2, 2 *(1 << j) + 1]);
        }
    }
    else if benchmark == "hist" {
        let k: usize = args[5].parse().unwrap();
        for j in first_t..=last_t {
            test_set.push([1<<j, 36*k*(1 << j), 96*k*(1 << j) - 1]);
        }
    }
    else if benchmark == "mlr" {
        //n=t/4 & $2n^2+111n-282$ &  $2n^2+111n-282$\\
        //n(3k+2k^2)+k^3+2k^2+3k constraints,  and ≤n(9/2 k+ 3/2 k^2)+k^3+2k^2+3k variables. 
        let k: usize = args[5].parse().unwrap();
        for j in first_t..=last_t {
            //test_set.push([4*(1<<j), 2*(1 << (j*2))+111*(1 << j) - 282, 2*(1 << (j*2))+111*(1 << j) - 283]);
            test_set.push([4*(1<<j), (1 << j)*(3*k + 2*k*k)+k*k*k + 2*k*k + 3*k, (1 << j)*(4*k + k*k)+k*k*k + 2*k*k + 3*k]);
        }
    }
    else {
        let fixedsize: usize = args[5].parse().unwrap();
        for j in first_t..=last_t {
            test_set.push([1<<j, (1 << fixedsize), (1 << fixedsize) - 1]);
        }
    }
/*
    //set of parameters to test
    let test_set = vec![
            [1 << 5, 1 << 10 ],  //t=32, |H|=4*t
            [1 << 6, 1 << 16],   //t=32, |H|=t^2
        // Fixed computation size, growing inputs
        /*
            [1 << 8, 1 << 20],   //t=256, |H|~1M
            [1 << 10, 1 << 20],   //t=1024, |H|~1M
            [1 << 12, 1 << 20],   //t=4096, |H|~1M
            [1 << 14, 1 << 20],   //t=16384, |H|~1M
            [1 << 16, 1 << 20],   //t~64K, |H|~1M
        */
        // Fixed computation medium size, growing inputs
        /*
            [1 << 7, 1 << 16],   //t=128, |H|~64K
            [1 << 8, 1 << 16],   //t=256, |H|~64K
            [1 << 9, 1 << 16],   //t=512, |H|~64K
            [1 << 10, 1 << 16],   //t=1024, |H|~64K
            [1 << 11, 1 << 16],   //t=2048, |H|~64K
            [1 << 12, 1 << 16],   //t=4096, |H|~64K
            [1 << 13, 1 << 16],   //t=8192, |H|~64K
            [1 << 14, 1 << 16],   //t=16384, |H|~64K
        */
        // Growing inputs, computation linear in input size
        /*    
            [1 << 7, 1 << 9],   //t=128, |H|=4*t
            //[1 << 8, 1 << 10],   //t=256, |H|=4*t
            [1 << 9, 1 << 11],   //t=512, |H|=4*t
            //[1 << 10, 1 << 12],   //t=1024, |H|=4*t
            [1 << 11, 1 << 13],   //t=2048, |H|=4*t
            //[1 << 12, 1 << 14],   //t=4096, |H|=4*t
            [1 << 13, 1 << 15],   //t=8192, |H|=4*t
            //[1 << 14, 1 << 16],   //t=16384, |H|=4*t
            [1 << 15, 1 << 17],   //t=16384, |H|=4*t
            //[1 << 16, 1 << 18],   //t~64K, |H|=4*t
        */
        // Growing inputs, computation quadratic in input size
        /*    
            [1 << 7, 1 << 14],   //t=128, |H|=t^2
            [1 << 8, 1 << 16],   //t=256, |H|=t^2
            [1 << 9, 1 << 18],   //t=512, |H|=t^2
            [1 << 10, 1 << 20],   //t=1024, |H|=t^2
        */
    ];
    */
    eprintln!("tests = {:?}", test_set);
    println!("HSNP '{}' benchmarks for {} with precomputed hash={} ({} repetitions) in microsecs:\nt, constraints, variables, signing-time, proving-time, ver-time",
        benchmark, stringify!(Bls12_381), precompute_hash, NUM_REPETITIONS);
    /*
    for test in test_set.iter() {
        let t = test[0];
        let num_cnstr = test[1];
        
        hsnp_benchmark!(bls, Fr, Bls12_381, t, num_cnstr);    
    }*/

    hsnp_benchmark!(bls, Fr, Bls12_381, precompute_hash, test_set);
    /*
    let t: usize = 1 << 5;
    let num_cnstr: usize = 4*t;
    
    hsnp_benchmark!(bls, Fr, Bls12_381, t, num_cnstr);
    */
}

#[derive(Copy)]
struct DummyCircuit<F: PrimeField> {
    pub a: Option<F>,
    pub b: Option<F>,
    pub num_variables: usize,
    pub num_constraints: usize,
}

impl<F: PrimeField> Clone for DummyCircuit<F> {
    fn clone(&self) -> Self {
        DummyCircuit {
            a: self.a.clone(),
            b: self.b.clone(),
            num_variables: self.num_variables.clone(),
            num_constraints: self.num_constraints.clone(),
        }
    }
}

impl<F: PrimeField> ConstraintSynthesizer<F> for DummyCircuit<F> {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        let a = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
        let b = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;
        let c = cs.new_input_variable(|| {
            let a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
            let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;

            Ok(a * b)
        })?;

        for _ in 0..(self.num_variables - 3) {
            let _ = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
        }

        for _ in 0..self.num_constraints - 1 {
            cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;
        }

        cs.enforce_constraint(lc!(), lc!(), lc!())?;

        Ok(())
    }
}