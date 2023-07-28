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

// don't use From trait on slices, unless you really mean slices
impl From<&Vec<Module>> for HMM {
    fn from(modules: &Vec<Module>) -> Self {
        let mut n_nodes = 0;
        for m in modules {
            n_nodes += match m {
                Module::Sequence(x) => { x.len() },
                Module::Repeat(x) => { x.0.len() },
            }
        }
        // starting background + flanks and motif + ending background + inserts
        n_nodes = 1 + n_nodes + 1 + (n_nodes-1);

        let mut model = HMM {
            initial: Array::zeros(n_nodes),
            transition: Array::zeros((n_nodes, n_nodes)),
            emission: Array::zeros((5, n_nodes))
        };

        model
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

    #[test]
    fn construct_hmm() {
        let modules: Vec<Module> = vec![
            (&b"TCT"[..]).into(),   // inside parenthesis creates &[u8] instead of &[u8;N]
            (&b"TAC"[..], 8).into(),
            (&b"GTC"[..], 5).into(),
            (&b"AAA"[..]).into()
        ];

        let model = HMM::from(&modules);
    }
}

