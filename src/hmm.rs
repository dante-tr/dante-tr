use std::result;

use ndarray::{Array, stack, ArrayView1, ArrayBase};
use ndarray;
use ndarray::s;
use ndarray_npy::read_npy;
use ndarray::Axis;
// https://docs.rs/ndarray/latest/ndarray/doc/ndarray_for_numpy_users/index.html

// maybe I should use this?
// https://docs.rs/hmmm/latest/hmmm/struct.HMM.html

const QUALITY_START: u8 = 33;
const N_QUAL: usize = 94;
const N_NUCL: usize = 5;
const NUCLEOTIDE_INDEX: [usize; 256] = {
    let mut map = [5; 256];
    map[b'A' as usize] = 0;
    map[b'C' as usize] = 1;
    map[b'G' as usize] = 2;
    map[b'T' as usize] = 3;
    map[b'N' as usize] = 4;
    map
};
const BASE_N_PROB: f32 = 0.001;
const P_INS: f32 = 1e-4;
const P_DEL: f32 = 1e-4;
const P_SNP: f32 = 0.0005;
const FREQ: f32 = 0.001;

#[derive(Default)]
pub struct HMM {
    pub initial: ndarray::Array1<f32>,
    pub transition: ndarray::Array2<f32>,
    pub emission: ndarray::Array3<f32>,
}

enum State {
    StartBackground,
    Sequence,
    SequenceEnd,
    Motif,
    MotifEnd,
    EndBackground,
    Insert,
}

// don't use From trait on slices, unless you really mean slices
impl From<&Vec<Module>> for HMM {
    fn from(modules: &Vec<Module>) -> Self {
        let motif_frequency = 0.001;
        let p_ins = 1e-4;
        let p_del = 1e-4;

        let seq = b"XACTGTGCAGXXXXXXXXXX".to_vec();
        let mut states: Vec<State> = Vec::new();
        states.push(State::StartBackground);
        for m in modules {
            match m {
                Module::Sequence(x) => {
                    for _ in 0..x.len()-1 {
                        states.push(State::Sequence);
                    }
                    states.push(State::SequenceEnd);
                },
                Module::Repeat(x) => {
                    for _ in 0..x.0.len()-1 {
                        states.push(State::Motif);
                    }
                    states.push(State::MotifEnd);
                },
            }
        }
        states.push(State::EndBackground);
        // insert states are only between 2 non-background states (e.g. BSISISB)
        for _ in 1..states.len()-2 {
            states.push(State::Insert);
        }

        let initial = initial_probabilities(&states, motif_frequency);
        let transition = transition_probabilities(
            &states, motif_frequency, p_del, p_ins, modules
        );
        let emission = emission_probabilities(&states, &seq);

        HMM { initial, transition, emission }
    }
}

fn initial_probabilities(states: &Vec<State>, init: f32) -> ndarray::Array1<f32> {
    let mut probabilities = Array::zeros(states.len());
    let mut n = 0;
    for (i, state) in states.iter().enumerate() {
        match state {
            State::Sequence | State::SequenceEnd | State::Motif | State::MotifEnd 
                => {
                    probabilities[[i]] = init;
                    n += 1; 
                },
            _ => {},
        }
    }
    probabilities[[0]] = 1.0 - (init * n as f32);
    assert!(probabilities[[0]] >= 0.0);
    return probabilities;
}

fn transition_probabilities(
    states: &Vec<State>,
    mfreq: f32,
    p_del: f32,
    p_ins: f32,
    modules: &Vec<Module>
) -> ndarray::Array2<f32> {
    let mut p = Array::zeros((states.len(), states.len()));

    let bg_start = 0;
    let module_starts = vec![1, 4, 7];
    let module_end = vec![3, 6, 9];
    let bg_end = 10;
    let k = 2;

    // connect first state
    // connect background start
    for i in module_starts.iter() { p[[bg_start, *i]] = mfreq; }
    p[[bg_start, bg_start + k]] = mfreq * p_del;
    p[[bg_start, bg_start]] = 1.0 - mfreq * module_starts.len() as f32 - mfreq * p_del;

    // connect deletions
    //  connect linear
    for i in 1..bg_end-k+1 {
        p[[i, i+k]] = p_del;
    }

    for i in 0..module_starts.len() {
        if matches!(states[module_starts[i]], State::Sequence) { continue; } 
        let l = module_end[i]-module_starts[i]+1;
        for j in 0..l {
            let del_start = module_starts[i] + j;
            let del_end = module_starts[i] + (j + k) % l;
            p[[del_start, del_end]] = p_del;
        }
    }

    // connect insertions
    for i in 1..bg_end-2+1 {
        let ins = bg_end + i;
        p[[i, ins]] = p_ins;
        p[[ins, ins]] = p_ins;
        if matches!(states[i], State::MotifEnd) {
            let (m_start, rep) = (4, 5);
            p[[ins, m_start]] = 1.0 - 1.0/rep as f32;
        }
        let used = p.slice(s![ins, ..]).sum();
        p[[ins, i+1]] = 1.0 - used;
    }
    // connect cycles
    p[[module_end[1], module_starts[1]]] = 1.0 - 1.0/5 as f32;

    // connect modules
    for (i, end) in module_end.iter().enumerate() {
        let used = p.slice(s![*end, ..]).sum();
        let prob = (1.0 - used)/((module_starts[i+1..].len() + 1) as f32);
        for start in module_starts[i+1..].iter() {
            p[[*end, *start]] = prob;
        }
        p[[*end, bg_end]] = prob;
    }

    // connect next
    for i in 1..bg_end-1 {
        let used = p.slice(s![i, ..]).sum();
        if used < 1.0 {
            p[[i, i+1]] = 1.0 - used;
        }
    }

    // connect last state
    p[[bg_end, bg_end]] = 1.0;

    return p;
}

fn emission_probabilities(states: &Vec<State>, letters: &Vec<u8>) -> ndarray::Array3<f32> {
    let mut result = Array::zeros((N_NUCL, N_QUAL, states.len()));

    for i in 0..N_NUCL {
        for j in 0..N_QUAL {
            for k in 0..states.len() {
                match states[k] {
                    State::StartBackground | State::EndBackground | State::Insert
                    => {
                        if i == NUCLEOTIDE_INDEX[b'N' as usize] {
                            result[[i, j, k]] = BASE_N_PROB; 
                        } else { result[[i, j, k]] = (1.0 - BASE_N_PROB) / 4.0; }
                    }
                    State::Sequence | State::SequenceEnd | State::Motif | State::MotifEnd
                    => {
                        let p_correct = 1.0 - P_SNP - 10.0f32.powf(-(j as f32)/10.0) - BASE_N_PROB;
                        let p_correct = p_correct.max((1.0 - BASE_N_PROB) / 4.0);
                        if i == NUCLEOTIDE_INDEX[letters[k] as usize] {
                            result[[i, j, k]] = p_correct;
                        } else if i == NUCLEOTIDE_INDEX[b'N' as usize] {
                            result[[i, j, k]] = BASE_N_PROB;
                        } else {
                            result[[i, j, k]] = (1.0 - p_correct - BASE_N_PROB) / 3.0;
                        }
                    }
                }
            }
        }
    }

    let result = result.map(|&x| x.ln());
    return result;
}

fn argmax(a: ArrayView1<f32>) -> usize {
    let mut max_pos = 0;
    let mut max_val = f32::NEG_INFINITY;
    for (p, &v) in a.iter().enumerate() {
        if v > max_val { max_pos = p; max_val = v; }
    }
    return max_pos;
}

impl HMM {
    /// Returns the most likely path and its log-likelihood of (seq, qual) given HMM.
    /// Assumes HMM with more than 0 states and probabilities in log space.
    pub fn log_predict(&self, seq: &[u8], qual: &[u8]) -> (f32, Vec<usize>) {
        let seq: Vec<_> = seq.iter().map(|x| NUCLEOTIDE_INDEX[*x as usize]).collect();
        let qual: Vec<_> = qual.iter().map(|x| (x - QUALITY_START) as usize).collect();

        let n_states = self.initial.len();
        let mut backptr = Array::ones((seq.len(), n_states)) * usize::MAX;

        let mut trellis: ndarray::Array1<f32> = &self.initial
            + &self.emission.slice(s![seq[0], qual[0], ..]);

        for i in 1..seq.len() {
            // for each pair of states transit
            let x = stack(Axis(1), &vec![trellis.view(); n_states]).unwrap() // safe
                + &self.transition;

            // for each state find incoming connection
            let incoming = x.map_axis(Axis(0), argmax);
            backptr.slice_mut(s![i, ..]).assign(&incoming);

            // for each state calculate probability of emitting (seq[i], qual[i])
            // next line mirrors `trellis = x.map_axis(Axis(0), max)` but in O(n)
            for j in 0..n_states { trellis[j] = x[[incoming[j], j]]; }
            trellis += &self.emission.slice(s![seq[i], qual[i], ..]);
        }

        let best_end = argmax(trellis.view());
        let likelihood = trellis[best_end];

        let mut path = vec![n_states; seq.len()];
        let last = path.last_mut().unwrap(); // safe
        *last = best_end;
        for i in (0..seq.len()-1).rev() {
            path[i] = backptr[[i+1, path[i+1]]];
        }

        (likelihood, path)
    }
}

enum Module {
    Sequence(Vec<u8>),
    Repeat((Vec<u8>, usize)),
}

impl From<&[u8]> for Module {
    fn from(s: &[u8]) -> Self {
        Module::Sequence(s.to_vec())
    }
}

impl From<(&[u8], usize)> for Module {
    fn from(tuple: (&[u8], usize)) -> Self {
        Module::Repeat((tuple.0.to_vec(), tuple.1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::approx_eq;
    use ndarray::ArrayView;
    use ndarray::Dim;
    use ndarray_npy::read_npy;
    use ndarray::Array2;
    use ndarray::Array3;

    #[test]
    fn construct_hmm() {
        let modules: Vec<Module> = vec![
            (&b"TCT"[..]).into(),   // inside parenthesis creates &[u8] instead of &[u8;N]
            // (&b"TAC"[..], 8).into(),
            (&b"GTC"[..], 5).into(),
            (&b"AAA"[..]).into()
        ];

        let model = HMM::from(&modules);
        let expected: Array2<f32> = read_npy("data/log_trans_f32.npy").unwrap();

        for i in 0..expected.shape()[0] { for j in 0..expected.shape()[1] {
            assert!(
                approx_eq!(f32, expected[[i, j]], model.transition[[i, j]].ln(), (1e-3, 2)),
                "for i={i} and j={j}: Expected {}, got {}.", 
                expected[[i, j]], model.transition[[i, j]]
            );
        }}

        // println!("{:#?}", model.initial);
        // println!("{:#?}", model.transition);
        // println!("{:#?}", model.emission);
    }

    #[test]
    fn prediction_works() {
        let initial = read_npy("data/log_init_f32.npy").unwrap();
        let transition = read_npy("data/log_trans_f32.npy").unwrap();
        let emission = read_npy("data/log_emit_f32.npy").unwrap();

        let model = HMM { initial, transition, emission };
        println!("Initial strides: {:?}", model.initial.strides());
        println!("Initial shape: {:?}", model.initial.shape());
        println!("Transition strides: {:?}", model.transition.strides());
        println!("Transition shape: {:?}", model.transition.shape());
        println!("Emission strides: {:?}", model.emission.strides());
        println!("Emission shape: {:?}", model.emission.shape());

        let seq =  b"AATCTGTCGTCGTCGTCAGTCGTCAAATT".to_vec();
        let qual = b":F::FF:,F,FFFFFFF,FF,FFF:F,FF".to_vec();
        let (likelihood, path) = model.log_predict(&seq, &qual);
        println!("{likelihood}: {path:?}");

        assert!(approx_eq!(f32, likelihood, 7.106122e-13_f32.ln(), (1e-3, 2)));
        assert!(path == vec![
            0, 0, 1, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4,
            5, 6, 16, 4, 5, 6, 4, 5, 6, 7, 8, 9, 10, 10
        ]);
    }

    #[test]
    fn emissions_are_correct() {
        let states = vec![
            State::StartBackground,
            State::Sequence, State::Sequence, State::SequenceEnd,
            State::Motif, State::Motif, State::MotifEnd,
            State::Sequence, State::Sequence, State::SequenceEnd,
            State::EndBackground,
            State::Insert, State::Insert, State::Insert, State::Insert,
            State::Insert, State::Insert, State::Insert, State::Insert
        ];
        let seq = b"_TCTGTCAAA_IIIIIIII".to_vec();
        let obtained = emission_probabilities(&states, &seq);
        let expected: Array3<f32> = read_npy("data/log_emit_f32.npy").unwrap();
        assert_eq_ndarray3(expected.view(), obtained.view(), (1e-3, 2));
    }

    // how to make this generic over dimensions?
    fn assert_eq_ndarray3(
        a1: ArrayView<f32, Dim<[usize; 3]>>,
        a2: ArrayView<f32, Dim<[usize; 3]>>, 
        acc: (f32, i32)
    ) {
        let shp = a1.shape();
        for i in 0..shp[0] {
            for j in 0..shp[1] {
                for k in 0..shp[2] {
                    assert!(
                        approx_eq!(f32, a1[[i, j, k]], a2[[i, j, k]], acc),
                        "for i={i}, j={j}, k={k}: Expected {}, got {}.",
                        a1[[i, j, k]].exp(), a2[[i, j, k]].exp()
                    )
                }
            }
        }
    }
}

