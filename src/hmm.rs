use ndarray::Array;
use ndarray;

// https://docs.rs/ndarray/latest/ndarray/
// maybe I should use this?
// https://docs.rs/hmmm/latest/hmmm/struct.HMM.html

#[derive(Default)]
pub struct HMM {
    pub initial: ndarray::Array1<f32>,
    pub transition: ndarray::Array2<f32>,
    pub emission: ndarray::Array2<f32>,
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
        // let transition = transition_probabilities(&mut model.transition);
        // let emission = emission_probabilities(&mut model.emission);

        HMM {
            initial,
            transition: Array::zeros((states.len(), states.len())),
            emission: Array::zeros((5, states.len()))
        }
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

fn transition_probabilities(probs: &mut ndarray::Array2<f32>) {}
fn emission_probabilities(probs: &mut ndarray::Array2<f32>) {}

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

    #[test]
    fn construct_hmm() {
        let modules: Vec<Module> = vec![
            (&b"TCT"[..]).into(),   // inside parenthesis creates &[u8] instead of &[u8;N]
            (&b"TAC"[..], 8).into(),
            (&b"GTC"[..], 5).into(),
            (&b"AAA"[..]).into()
        ];

        let model = HMM::from(&modules);
        println!("{}", model.initial);
    }
}

