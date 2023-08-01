use ndarray::{Array, stack};
use ndarray;
use ndarray::s;
use ndarray_npy::read_npy;
use ndarray::Axis;
// https://docs.rs/ndarray/latest/ndarray/doc/ndarray_for_numpy_users/index.html

// maybe I should use this?
// https://docs.rs/hmmm/latest/hmmm/struct.HMM.html

const QUALITY_START: u8 = 33;
const NUCLEOTIDE_INDEX: [u8; 256] = {
    let mut map = [5; 256];
    map[b'A' as usize] = 0;
    map[b'C' as usize] = 1;
    map[b'G' as usize] = 2;
    map[b'T' as usize] = 3;
    map[b'N' as usize] = 4;
    map
};
const BASE_N_PROB: f32 = 0.001;

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
        let emission = emission_probabilities();

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

fn emission_probabilities() -> ndarray::Array3<f32> {
    let emissions: ndarray::Array3<f32> = read_npy("data/emissions.npy").unwrap();
    return emissions;
}

impl HMM {
    fn predict(&self, seq: &[u8], qual: &[u8]) -> (f32, Vec<usize>) {
        let mut backptr = Array::ones((self.initial.len(), seq.len())) * -1.0;
        let tmp: f32 = backptr[[0, 0]];
        println!("{:?} {}", backptr.shape(), tmp);

        println!("{:?}", self.emission.shape());
        let mut trellis: ndarray::Array1<f32> = &self.initial
            + &self.emission.slice(s![seq[0] as usize, qual[0] as usize, ..]);
        println!("{:?}", trellis.shape());

        for i in 1..seq.len() {
            // for each pair of states transit
            let tmp = stack(Axis(0), &vec![trellis.view(); self.transition.shape()[0]]).unwrap();
            let tmp = &tmp + &self.transition;

            // for each state find incoming connection
            let tmp2 = tmp.sum_axis(Axis(1)); // should be argmax
            backptr.slice_mut(s![.., i]).assign(&tmp2);

            // for each state calculate probability of emitting (seq[i], qual[i])
            // max + emissions
            trellis = &tmp2 + &self.emission.slice(s![seq[i] as usize, qual[i] as usize, ..]);
        }

        //     temp_mat = (
        //         np.tile(trellis[t - 1, :], (self.num_states, 1)) + self.trans_prob
        //     )
        //     backpt[t, :] = np.argmax(temp_mat, axis=1)
        //     trellis[t, :] = (
        //         temp_mat[np.arange(self.num_states), backpt[t, :]]
        //         + self.obs_prob[obs[t], quality[t]]
        //     )

        let likelihood = trellis.sum();

        let mut path = vec![self.initial.shape()[0]+1; seq.len()];
        let last = path.last_mut().unwrap();
        *last = trellis.map_axis(Axis(0), |_| 0).sum(); // should be argmax
        for i in (0..seq.len()-1).rev() {
            path[i] = 0; // backptr[[i, path[i+1]]];
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
    use ndarray_npy::read_npy;
    use ndarray::Array2;

    #[test]
    fn construct_hmm() {
        let modules: Vec<Module> = vec![
            (&b"TCT"[..]).into(),   // inside parenthesis creates &[u8] instead of &[u8;N]
            // (&b"TAC"[..], 8).into(),
            (&b"GTC"[..], 5).into(),
            (&b"AAA"[..]).into()
        ];

        let model = HMM::from(&modules);
        let expected: Array2<f32> = read_npy("data/transitions_f32.npy").unwrap();

        for i in 0..expected.shape()[0] { for j in 0..expected.shape()[1] {
            assert!(
                approx_eq!(f32, expected[[i, j]], model.transition[[i, j]], (1e-3, 2)),
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
        let modules: Vec<Module> = vec![
            (&b"TCT"[..]).into(),
            (&b"GTC"[..], 5).into(),
            (&b"AAA"[..]).into(),
        ];

        let model = HMM::from(&modules);

        let seq =  b"AATCTGTCGTCGTCGTCAGTCGTCAAATT".to_vec();
        let qual = b":F::FF:,F,FFFFFFF,FF,FFF:F,FF".to_vec();
        let seq: Vec<_> = seq.iter().map(|x| NUCLEOTIDE_INDEX[*x as usize]).collect();
        let qual: Vec<_> = qual.iter().map(|x| x - QUALITY_START).collect();
        println!("{seq:?}");
        println!("{qual:?}");
        let (likelihood, path) = model.predict(&seq, &qual);
        println!("{likelihood} {path:?}");
        // 7.106122e-13
        // [0, 0, 1, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6, 16, 4, 5, 6, 4, 5, 6, 7, 8, 9, 10, 10]
    }
}

