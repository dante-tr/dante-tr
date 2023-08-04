use ndarray::{self, s, Axis, Array, stack, ArrayView1};
// https://docs.rs/ndarray/latest/ndarray/doc/ndarray_for_numpy_users/index.html

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
const P_BASE_N: f32 = 0.001;
const P_INS: f32 = 1e-4;
const P_DEL: f32 = 1e-4;
const P_SNP: f32 = 0.0005;
const FREQ: f32 = 0.001;
const DEL: usize = 2; // represents deletion of 1 nucleotide

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

struct MDesc {
    start: usize,
    len: usize,
    rep: Option<usize>
}

enum State {
    Start,
    Seq{nucl: u8},
    Motif{nucl: u8},
    End,
    Ins,
}

#[derive(Default)]
pub struct HMM {
    pub initial: ndarray::Array1<f32>,
    pub transition: ndarray::Array2<f32>,
    pub emission: ndarray::Array3<f32>,
}

impl From<&Vec<Module>> for HMM {
    fn from(modules: &Vec<Module>) -> Self {
        let states = get_states(modules);
        let description = get_description(modules);

        let initial = initial_probabilities(&states);
        let transition = transition_probabilities(&states, &description);
        let emission = emission_probabilities(&states);

        HMM { initial, transition, emission }
    }
}

fn get_states(modules: &Vec<Module>) -> Vec<State> {
    let mut states: Vec<State> = Vec::new();
    states.push(State::Start);
    for module in modules {
        match module {
            Module::Sequence(x) => {
                for &c in x { states.push(State::Seq{nucl: c}); }
            },
            Module::Repeat((x, _)) => {
                for &c in x { states.push(State::Motif{nucl: c}); }
            }
        }
    }
    states.push(State::End);
    // insert states are only between 2 non-background states (e.g. BSISISB)
    for _ in 1..states.len()-2 { states.push(State::Ins); }
    return states;
}

fn get_description(modules: &Vec<Module>) -> Vec<MDesc> {
    let mut start = 1;
    let mut description = Vec::new();
    for module in modules {
        match module {
            Module::Sequence(s) => {
                description.push(MDesc{ start, len: s.len(), rep: None });
                start += s.len();
            },
            Module::Repeat((s, rep)) => {
                description.push(MDesc{ start, len: s.len(), rep: Some(*rep) });
                start += s.len();
            }
        }
    }
    return description;
}

fn initial_probabilities(states: &Vec<State>) -> ndarray::Array1<f32> {
    let mut p = Array::zeros(states.len());
    let mut n = 0;
    for (i, state) in states.iter().enumerate() {
        match state {
            State::Seq{..} | State::Motif{..} => {
                p[[i]] = FREQ;
                n += 1; 
            },
            _ => {},
        }
    }
    p[[0]] = 1.0 - (FREQ * n as f32);
    assert!(p[[0]] >= 0.0);
    return p;
}

fn transition_probabilities(states: &Vec<State>, desc: &Vec<MDesc>) -> ndarray::Array2<f32> {
    let mut p = Array::zeros((states.len(), states.len()));

    let bg_start = 0;
    let mut bg_end = 0;
    while !matches!(states[bg_end], State::End) { bg_end += 1; }

    // create intramodule connections
    for m in desc.iter() {
        if matches!(states[m.start], State::Seq{..}) { continue; }
        let m_end = m.start + m.len - 1;
        let m_rep = m.rep.unwrap() as f32; // safe due to previous if

        // add cycle
        p[[m_end, m.start]] = 1.0 - 1.0/m_rep;

        // add insertion between repetitions
        let ins = bg_end + m_end;
        p[[ins, m.start]] = 1.0 - 1.0/m_rep - P_INS;

        // add deletions between repetitions
        for i in 0..m.len {
            let del_start = m.start + i;
            let del_end = m.start + (i + DEL) % m.len;
            p[[del_start, del_end]] = P_DEL;
        }
    }

    // connect simple insertions
    for i in 1..=bg_end-2 {
        let ins = bg_end + i;
        p[[i, ins]] = P_INS;
        p[[ins, ins]] = P_INS;
        p[[ins, i+1]] = 1.0 - p.slice(s![ins, ..]).sum();
    }

    // connect simple deletions
    p[[bg_start, bg_start + DEL]] = P_DEL * FREQ;
    for i in 1..=bg_end-DEL {
        p[[i, i+DEL]] = P_DEL;
    }

    // allow module skipping
    for m in desc.iter() {
        p[[bg_start, m.start]] = FREQ;
    }

    for (i, m) in desc.iter().enumerate() {
        let m_end = m.start + m.len - 1;
        // number of remaining modules to the right
        let r_mod = (desc.len() - i) as f32;
        let r_prob = 1.0 - p.slice(s![m_end, ..]).sum();
        // +2, because we jump over the next module, as it will be connected later
        if let Some(x) = desc.get(i+2..) {
            for destination in x {
                p[[m_end, destination.start]] = r_prob / r_mod;
            }
        }
        p[[m_end, bg_end]] = r_prob / r_mod;
    }

    // loop in start
    p[[bg_start, bg_start]] = 1.0 - p.slice(s![bg_start, ..]).sum();

    // connect to the next state
    for i in 1..bg_end-1 {
        p[[i, i+1]] = 1.0 - p.slice(s![i, ..]).sum();
    }

    // loop in end
    p[[bg_end, bg_end]] = 1.0;

    return p;
}

fn emission_probabilities(states: &Vec<State>) -> ndarray::Array3<f32> {
    let mut p = Array::zeros((N_NUCL, N_QUAL, states.len()));

    for i in 0..N_NUCL { for j in 0..N_QUAL { for k in 0..states.len() {
        match states[k] {
            State::Start | State::End | State::Ins => {
                if i == NUCLEOTIDE_INDEX[b'N' as usize] {
                    p[[i, j, k]] = P_BASE_N; 
                } else {
                    p[[i, j, k]] = (1.0 - P_BASE_N) / 4.0;
                }
            }
            State::Seq{nucl: c} | State::Motif{nucl: c} => {
                let p_correct = 1.0 - P_SNP - 10.0f32.powf(-(j as f32)/10.0) - P_BASE_N;
                let p_correct = f32::max(p_correct, (1.0 - P_BASE_N) / 4.0);
                if i == NUCLEOTIDE_INDEX[c as usize] {
                    p[[i, j, k]] = p_correct;
                } else if i == NUCLEOTIDE_INDEX[b'N' as usize] {
                    p[[i, j, k]] = P_BASE_N;
                } else {
                    p[[i, j, k]] = (1.0 - p_correct - P_BASE_N) / 3.0;
                }
            }
        }
    }}}

    let result = p.map(|&x| x.ln());
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
            (&b"GTC"[..], 5).into(),
            (&b"AAA"[..]).into()
        ];

        let model = HMM::from(&modules);
        let expected: Array2<f32> = read_npy("data/log_trans_f32.npy").unwrap();
        let expected = expected.map(|&x| x.exp());

        let diff = find_diff_ndarray2(expected.view(), model.transition.view(), (1e-4, 2));
        if let Some((i, j)) = diff {
            println!("{} {} {} {}", i, j, expected[[i, j]], model.transition[[i, j]]);
            println!("{:#034b}", expected[[i, j]].to_bits());
            println!("{:#034b}", model.transition[[i, j]].to_bits());
        }
        assert!(diff.is_none());

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

//     #[test]
//     fn initial_are_correct() {
//         let states = vec![
//             State::StartBackground,
//             State::Sequence, State::Sequence, State::SequenceEnd,
//             State::Motif, State::Motif, State::MotifEnd,
//             State::Sequence, State::Sequence, State::SequenceEnd,
//             State::EndBackground,
//             State::Insert, State::Insert, State::Insert, State::Insert,
//             State::Insert, State::Insert, State::Insert, State::Insert
//         ];
//         let obtained = initial_probabilities(&states);
//         let expected = read_npy("data/log_init_f32.npy").unwrap();
//     }

    #[test]
    fn emissions_are_correct() {
        let states = vec![
            State::Start,
            State::Seq{nucl: b'T'}, State::Seq{nucl: b'C'}, State::Seq{nucl: b'T'},
            State::Motif{nucl: b'G'}, State::Motif{nucl: b'T'}, State::Motif{nucl: b'C'},
            State::Seq{nucl: b'A'}, State::Seq{nucl: b'A'}, State::Seq{nucl: b'A'},
            State::End,
            State::Ins, State::Ins, State::Ins, State::Ins,
            State::Ins, State::Ins, State::Ins, State::Ins
        ];
        let obtained = emission_probabilities(&states);
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
            let diff = find_diff_ndarray2(
                a1.index_axis(Axis(0), i),
                a2.index_axis(Axis(0), i),
                acc
            );
            if let Some((j, k)) = diff {
                println!(
                    "for i={i}, j={j}, k={k}: Expected {}, got {}.",
                    a1[[i, j, k]], a2[[i, j, k]]
                );
                panic!();
            }
        }
    }

    fn find_diff_ndarray2(
        a1: ArrayView<f32, Dim<[usize; 2]>>,
        a2: ArrayView<f32, Dim<[usize; 2]>>,
        acc: (f32, i32)
    ) -> Option<(usize, usize)> {
        let shp = a1.shape();
        for i in 0..shp[0] {
            for j in 0..shp[1] {
                if !approx_eq!(f32, a1[[i, j]], a2[[i, j]], acc) {
                    return Some((i, j));
                }
            }
        }
        return None;
    }

    // fn find_diff<D>(a1: ArrayView<f32, D>, a2: ArrayView<f32, D>, acc: (f32, i32)) -> Option<Vec<usize>>
    // where
    //     D: Dimension + RemoveAxis,
    //     [usize; 1]: NdIndex<D>
    // {
    //     let shp = a1.shape();
    //     if shp.len() == 1 {
    //         for i in 0..shp[0] {
    //             if !approx_eq!(f32, a1[[i]], a2[[i]], acc) {
    //                 return Some(vec![i]);
    //             }
    //         }
    //         return None;
    //     } else {
    //         for i in 0..shp[0] {
    //             let tmp1 = a1.index_axis(Axis(0), i).into_dimensionality().unwrap();
    //             let tmp2 = a2.index_axis(Axis(0), i).into_dimensionality().unwrap();
    //             let x = find_diff(
    //                 tmp1,
    //                 tmp2,
    //                 acc
    //             );
    //             if let Some(mut x) = x {
    //                 let mut v = vec![i];
    //                 x.reverse();
    //                 v.extend_from_slice(&x);
    //                 return Some(v);
    //             }
    //         }
    //         return None;
    //     }
    // }

    // #[test]
    // fn test_find_diff() {
    //     let a1 = array![[1.,2.,3.], [4.,5.,6.]];
    //     let a2 = array![[1.,2.,3.], [4.,0.,6.]];

    //     // let x = find_diff(a1.view(), a2.view(), (1e-3, 2));
    // }
}

