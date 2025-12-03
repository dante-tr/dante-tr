use std::{collections::{HashMap, HashSet}, ops::Range};

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
const NUCLEOTIDE_BACKDEX: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];
const P_BASE_N: f32 = 0.001;
const P_INS: f32 = 1e-4;
const P_DEL: f32 = 1e-4;
const P_SNP: f32 = 0.0005;
const FREQ: f32 = 0.001;
const DEL: usize = 2; // represents deletion of 1 nucleotide

const ASCII_ZERO: u8 = 48;

pub enum Module {
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

#[derive(Debug)]
struct MDesc {
    start: usize,
    len: usize,
    rep: Option<usize>
}

#[derive(Debug)]
enum State {
    Start,
    Seq{c: u8, id: u8},
    Motif{c: u8, id: u8},
    End,
    Ins,
}

type TransitionSet = HashSet<(usize, usize)>;

#[derive(Default)]
pub struct Hmm {
    states: Vec<State>,
    state_to_mod: Vec<usize>,
    deletions: TransitionSet,
    module_changes: TransitionSet,
    unit_changes: TransitionSet,
    initial: ndarray::Array1<f32>,
    transition: ndarray::Array2<f32>,
    emission: ndarray::Array3<f32>,
}

impl From<&Vec<Module>> for Hmm {
    fn from(modules: &Vec<Module>) -> Self {
        let states = get_states(modules);
        let description = get_description(modules);
        let state_to_mod = create_state_to_mod(&states, &description);

        let initial = initial_probabilities(&states);
        let (transition, tsets) = transition_probabilities(&states, &description);
        let emission = emission_probabilities(&states);

        let (deletions, module_changes, unit_changes) = tsets;

        Hmm { states, state_to_mod, deletions, module_changes, unit_changes, initial, transition, emission }
    }
}

#[test]
fn test_construction_of_transition_sets() {
    let modules: Vec<Module> = vec![
        (&b"TTTT"[..]).into(),
        (&b"GCG"[..], 5).into(),
        (&b"TTTT"[..]).into()
    ];

    let states = get_states(&modules);
    let description = get_description(&modules);

    let (_, tsets) = transition_probabilities(&states, &description);

    let (trans_deletion, trans_modchange, trans_unitchange) = tsets;
    let exp_trans_deletion = HashSet::from([
        (0, 2), (1, 3), (2, 4), (3, 5), (4, 6), (5, 7), (6, 5), (6, 8), (7, 6),
        (7, 9), (8, 10), (9, 11), (10, 12)
    ]);
    assert!(exp_trans_deletion == trans_deletion);

    let exp_trans_modchange = HashSet::from([
        (0, 1), (0, 2), (0, 5), (0, 8), (3, 5), (4, 5), (4, 6), (4, 8), (4, 12),
        (6, 8), (7, 8), (7, 9), (7, 12), (10, 12), (11, 12), (16, 5), (19, 8)
    ]);
    assert!(exp_trans_modchange == trans_modchange);

    let exp_trans_unitchange = HashSet::from([(6, 5), (7, 5), (7, 6), (19, 5)]);
    assert!(exp_trans_unitchange == trans_unitchange);
}

fn get_states(modules: &[Module]) -> Vec<State> {
    let mut states: Vec<State> = Vec::new();
    states.push(State::Start);
    for (i, module) in modules.iter().enumerate() {
        let id = i as u8;
        match module {
            Module::Sequence(x) => {
                for &c in x { states.push(State::Seq{c, id}); }
            },
            Module::Repeat((x, _)) => {
                for &c in x { states.push(State::Motif{c, id}); }
            }
        }
    }
    states.push(State::End);
    // insert states are only between 2 non-background states (e.g. BSISISB)
    for _ in 1..states.len()-2 { states.push(State::Ins); }
    return states;
}

fn get_description(modules: &[Module]) -> Vec<MDesc> {
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

fn initial_probabilities(states: &[State]) -> ndarray::Array1<f32> {
    let mut p = Array::zeros(states.len());
    for (i, state) in states.iter().enumerate() {
        if matches!(state, State::Seq{..} | State::Motif{..}) { p[[i]] = FREQ; }
    }
    p[[0]] = 1.0 - p.sum();
    assert!(p[[0]] >= 0.0);
    return p;
}

fn create_state_to_mod(states: &[State], desc: &[MDesc]) -> Vec<usize> {
    let mut x = Vec::new();
    x.push(usize::MAX);  // starting in bg
    for (i, d) in desc.iter().enumerate() {
        x.extend_from_slice(&vec![i; d.len]);
    }

    let result = x.iter().cloned().cycle().take(states.len()).collect();
    return result;
}

/// Returns transition probability matrix and three transition sets (HashSet\<(usize, usize)\>) representing:
/// - transitions causing deletions
/// - transitions causing module change
/// - transitions causing unit change
fn transition_probabilities(states: &[State], desc: &[MDesc])
    -> (ndarray::Array2<f32>, (TransitionSet, TransitionSet, TransitionSet))
{
    let mut result = Array::zeros((states.len(), states.len()));
    let mut transition_delet = HashSet::new();  // transition causing deletion
    let mut transition_mchng = HashSet::new();  // transition causing module change
    let mut transition_uchng = HashSet::new();  // transition causing unit change

    let state_to_mod = create_state_to_mod(states, desc);

    let bg_start = 0;
    let mut bg_end = 0;
    while !matches!(states[bg_end], State::End) { bg_end += 1; }

    // create intramodule connections
    for m in desc.iter() {
        if matches!(states[m.start], State::Seq{..}) { continue; }
        let m_end = m.start + m.len - 1;
        let m_rep = m.rep.unwrap() as f32; // safe due to previous if

        // add cycle
        set_cell(&mut result, [m_end, m.start], 1.0 - 1.0/m_rep, "connection type 01");
        transition_uchng.insert((m_end, m.start));

        // add insertion between repetitions
        let ins = bg_end + m_end;
        if ins < states.len() {    // last -> bg does not have insertion
            let value = 1.0 - 1.0/m_rep - P_INS;  // connection type 02
            set_cell(&mut result, [ins, m.start], value, "connection type 02");
            transition_uchng.insert((ins, m.start));
        }

        // add deletions between repetitions
        for i in 0..m.len {
            let del_start = m.start + i;
            let del_end = m.start + (i + DEL) % m.len;
            if (del_start != m.start) & (del_end != m_end) {
                set_cell(&mut result, [del_start, del_end], P_DEL, "connection type 03");
                transition_uchng.insert((del_start, del_end));
                transition_delet.insert((del_start, del_end));
            }
        }
    }

    // connect simple insertions
    for i in 1..=bg_end-2 {
        let ins = bg_end + i;
        set_cell(&mut result, [i, ins], P_INS, "connection type 04");
        set_cell(&mut result, [ins, ins], P_INS, "connection type 05");

        let value = 1.0 - result.slice(s![ins, ..]).sum();
        set_cell(&mut result, [ins, i+1], value, "connection type 06");
        if state_to_mod[ins] != state_to_mod[i+1] { transition_mchng.insert((ins, i+1)); }
    }

    // connect simple deletions
    let cell = [bg_start, bg_start + DEL];
    set_cell(&mut result, cell, P_DEL * FREQ, "connection type 07");
    transition_delet.insert((bg_start, bg_start + DEL));
    transition_mchng.insert((bg_start, bg_start + DEL));

    for i in 1..=bg_end-DEL {
        set_cell(&mut result, [i, i + DEL], P_DEL, "connection type 08");
        transition_delet.insert((i, i + DEL));
        if state_to_mod[i] != state_to_mod[i+DEL] { transition_mchng.insert((i, i+DEL)); }
    }

    let m = &desc[0];
    set_cell(&mut result, [bg_start, m.start], FREQ, "start");
    transition_mchng.insert((bg_start, m.start));

    for (i, m) in desc.iter().enumerate().skip(desc.len() - 1) {
        let m_end = m.start + m.len - 1;
        // number of remaining modules to the right
        let r_mod = (desc.len() - i) as f32;
        let r_prob = 1.0 - result.slice(s![m_end, ..]).sum();
        set_cell(&mut result, [m_end, bg_end], r_prob / r_mod, "end");
        transition_mchng.insert((m_end, bg_end));
    }

    // loop in start
    let value = 1.0 - result.slice(s![bg_start, ..]).sum();
    set_cell(&mut result, [bg_start, bg_start], value, "connection type 12");

    // connect to the next state
    for i in 1..bg_end-1 {
        let value = 1.0 - result.slice(s![i, ..]).sum();  // connection type 13
        set_cell(&mut result, [i, i+1], value, "connection type 13");
        if state_to_mod[i] != state_to_mod[i+1] { transition_mchng.insert((i, i+1)); }
    }

    // loop in end
    set_cell(&mut result, [bg_end, bg_end], 1.0, "connection type 14");

    return (result, (transition_delet, transition_mchng, transition_uchng));
}

fn set_cell(matrix: &mut ndarray::Array2<f32>, cell: [usize; 2], value: f32, description: &str) {
    matrix[cell] = value;
    // println!("{:>2} -> {:>2}  {:>12}  {}", cell[0], cell[1], value, description);
    let _ = description;
}

fn emission_probabilities(states: &[State]) -> ndarray::Array3<f32> {
    let mut p = Array::zeros((N_NUCL, N_QUAL, states.len()));
    let n_pos = NUCLEOTIDE_INDEX[b'N' as usize];

    let iupac_to_nucls: HashMap<u8, HashSet<u8>> = HashMap::from([
        (b'A', HashSet::from([b'A'])),
        (b'C', HashSet::from([b'C'])),
        (b'G', HashSet::from([b'G'])),
        (b'T', HashSet::from([b'T'])),
        (b'M', HashSet::from([b'A', b'C'])),
        (b'R', HashSet::from([b'A', b'G'])),
        (b'W', HashSet::from([b'A', b'T'])),
        (b'S', HashSet::from([b'C', b'G'])),
        (b'Y', HashSet::from([b'C', b'T'])),
        (b'K', HashSet::from([b'G', b'T'])),
        (b'V', HashSet::from([b'A', b'C', b'G'])),
        (b'H', HashSet::from([b'A', b'C', b'T'])),
        (b'D', HashSet::from([b'A', b'G', b'T'])),
        (b'B', HashSet::from([b'C', b'G', b'T'])),
    ]); 

    for k in 0..states.len() { for i in 0..N_NUCL { for j in 0..N_QUAL {
        let cell = &mut p[[i, j, k]];

        use State as S;
        match states[k] {
            S::Start | S::End | S::Ins
            | S::Seq{c: b'N', ..} | S::Motif{c: b'N', ..}
            | S::Seq{c: b'X', ..} | S::Motif{c: b'X', ..} => {
                if i == n_pos { 
                    *cell = P_BASE_N;
                } else { 
                    *cell = (1.0 - P_BASE_N) / 4.0; 
                }
            }
            S::Seq{c, ..} | S::Motif{c, ..} => {
                let nucls = iupac_to_nucls.get(&c).expect("Unsupported.");
                let match_prob = p_c_eq_c_given_q(j, nucls.len());

                if i == n_pos { 
                    *cell = P_BASE_N;
                } else if nucls.contains(&NUCLEOTIDE_BACKDEX[i]) { 
                    *cell = match_prob;
                } else { 
                    *cell = (1.0 - P_BASE_N - match_prob) / (4 - nucls.len()) as f32}
            }
        }
    }}}

    return p;
}

/// TODO: write some explanation
/// Probability that letter c will be emitted. ???
fn p_c_eq_c_given_q(phred_score: usize, iupac_chars_number: usize) -> f32 {
    let p_base_calling_error = 10.0f32.powf(-(phred_score as f32)/10.0);
    let p_result =
        (1.0 - P_BASE_N - P_SNP - p_base_calling_error) / iupac_chars_number as f32;
    // let p_result = (1.0 - P_BASE_N - P_SNP - p_base_calling_error) / iupac_n;

    let p_random = (1.0 - P_BASE_N) / 4.0;
    return f32::max(p_result, p_random);
}

fn argmax(a: ArrayView1<f32>) -> usize {
    let mut max_pos = 0;
    let mut max_val = f32::NEG_INFINITY;
    for (p, &v) in a.iter().enumerate() {
        if v > max_val { max_pos = p; max_val = v; }
    }
    return max_pos;
}

impl Hmm {
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

    pub fn log(mut self) -> Self {
        self.initial = self.initial.map(|&x| x.ln());
        self.transition = self.transition.map(|&x| x.ln());
        self.emission = self.emission.map(|&x| x.ln());
        return self;
    }

    pub fn partition_to_units(&self, path: &[usize]) -> (Vec<Range<usize>>, Vec<usize>) {
        let mut partition = Vec::new();
        let mut mod_ids = Vec::new();
        let mut s = 0;

        for i in 1..path.len() {
            let mchng = self.module_changes.contains(&(path[i-1], path[i]));
            let uchng = self.unit_changes.contains(&(path[i-1], path[i]));
            if mchng || uchng {
                partition.push(s..i);
                mod_ids.push(self.state_to_mod[path[s]]);
                s = i;
            }
        }
        partition.push(s..path.len());
        mod_ids.push(self.state_to_mod[path[s]]);

        return (partition, mod_ids);
    }

    pub fn reconstruct_sequence(&self, path: &[usize]) -> Vec<u8> {
        let mut seq = Vec::with_capacity(path.len());
        for &i in path.iter() {
            match self.states[i] {
                State::Start | State::End => { seq.push(b'-'); },
                State::Ins => { seq.push(b'_'); }
                State::Seq{c, ..} | State::Motif{c, ..} => { seq.push(c); },
            }
        }
        return seq;
    }

    pub fn reconstruct_mod_ids(&self, path: &[usize]) -> Vec<u8> {
        let mut ids = Vec::with_capacity(path.len());
        for &i in path.iter() {
            match self.states[i] {
                State::Start | State::End => { ids.push(b'-'); },
                State::Ins => { ids.push(b'I'); },
                State::Seq {c:_, id} | State::Motif {c:_, id} => {
                    ids.push(ASCII_ZERO + id);
                }
            }
        }
        return ids;
    }

    pub fn realign(&self, path: &[usize], seq: &[u8]) -> (Vec<usize>, Vec<u8>) {
        let mut new_path = Vec::with_capacity(path.len());
        let mut new_seq = Vec::with_capacity(seq.len());

        for i in 0..path.len()-1 {
            new_path.push(path[i]);
            new_seq.push(seq[i]);
            // deletions should be maybe represented by HashMap
            // with keys (from, to) and values representing the nondeleted path
            if self.deletions.contains(&(path[i], path[i+1])) {
                new_path.push(path[i]+1);
                new_seq.push(b'_');
            }
        }
        new_path.push(*path.last().unwrap());
        new_seq.push(*seq.last().unwrap());
        return (new_path, new_seq);
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
    fn basic_test() {
        let modules: Vec<Module> = vec![
            (&b"TCT"[..]).into(),
            (&b"GTC"[..], 5).into(),
            (&b"AAA"[..]).into()
        ];

        let seq =  b"AATCTGTCGTCGTCGTCAGTCGTCAAATT".to_vec();
        let qual = b":F::FF:,F,FFFFFFF,FF,FFF:F,FF".to_vec();

        let model = Hmm::from(&modules).log();
        let (likelihood, annotation) = model.log_predict(&seq, &qual);

        assert!(approx_eq!(f32, likelihood, 7.106122e-13_f32.ln(), (1e-4, 2)));
        assert!(annotation == vec![
            0, 0, 1, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4,
            5, 6, 16, 4, 5, 6, 4, 5, 6, 7, 8, 9, 10, 10
        ]);

        let rebuilt = model.reconstruct_sequence(&annotation);
        let expected = b"--TCTGTCGTCGTCGTC_GTCGTCAAA--".to_vec();
        assert!(rebuilt == expected);
        let modules = model.reconstruct_mod_ids(&annotation);
        let expected = b"--000111111111111I111111222--";
        assert!(modules == expected);
    }

    #[test]
    fn reconstructed_output_has_the_same_length() {
        let modules: Vec<Module> = vec![
            (&b"ATACAAAAAAAAAAAAAAAA"[..]).into(),
            (&b"GAA"[..], 6).into(),
            (&b"AATAAAGAAAAGTTAGCCGG"[..]).into()
        ];
        let model = Hmm::from(&modules).log();

        let seq =  b"AAATAAAAAAAAAAAAAAAAAAAAGAAGAAGAAGAAGAAGAAGAAGAAGAAAATAAAGAAAAGTTAGCCGG".to_vec();
        let qual = b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF".to_vec();

        let (_, annotation) = model.log_predict(&seq, &qual);

        let (new_annot, read) = model.realign(&annotation, &seq); 
        let rref = model.reconstruct_sequence(&new_annot);
        let mods = model.reconstruct_mod_ids(&new_annot);

        assert_eq!(read.len(), rref.len());
        assert_eq!(rref, b"--ATACAAAAAAAAAAAAAAAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAAATAAAGAAAAGTTAGCCGG");
        assert_eq!(mods, b"--0000000000000000000011111111111111111111111111111122222222222222222222");
    }

    #[test]
    fn add_deletions_to_read() {
        let read =  b"AATTGTCGCGTCGTGTCGTCAAATT".to_vec();
        let annotation = vec![
            0, 0, 1, 3, 4, 5, 6, 4, 6, 4, 5, 6,
            4, 5, 4, 5, 6, 4, 5, 6, 7, 8, 9, 10, 10
        ];
        let expected = b"AAT_TGTCG_CGTCGT_GTCGTCAAATT".to_vec();

        let modules: Vec<Module> = vec![
            (&b"TCT"[..]).into(),
            (&b"GTC"[..], 5).into(),
            (&b"AAA"[..]).into()
        ];

        let model = Hmm::from(&modules).log();

        let (_, read) = model.realign(&annotation, &read);
        assert!(read == expected);
    }

    #[test]
    fn can_construct_from_single_module() {
        let modules = vec![
            (&b"ATTTT"[..], 30).into()
        ];
        let _model = Hmm::from(&modules).log();
    }

    #[test]
    fn predicts_single_letter_motifs() {
        let modules = vec![
            (&b"CTTGTTACTAAGCCTGATTT"[..]).into(),
            (&b"A"[..], 11).into(),
            (&b"TTACTTTCAGATGTCTGTCA"[..]).into(),    
        ];
        let model = Hmm::from(&modules).log();

        let sequence =
            b"AAGCCTGATTTAAAAAAAAAAAAAATTACTTTCAGATGT".to_vec();
        let quality = 
            b"FFFFFFFFFFFFFFF::FFF:FFFFFFF:FFFFFF:FFF".to_vec();
        let (_likelihood, annotation) = model.log_predict(&sequence, &quality);

        // assert!(
        //     approx_eq!(f32, likelihood, 4.019_535e-6_f32.ln(), (1e-3, 2)),
        //     "{likelihood} != {}", 4.019_535e-6_f32.ln()
        // );
        assert!(annotation == vec![
            10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
            22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35
        ]);
    }

    // #[ignore = "deprecated"]
    #[test]
    fn construct_hmm() {
        let modules: Vec<Module> = vec![
            (&b"TCT"[..]).into(),   // inside parenthesis creates &[u8] instead of &[u8;N]
            (&b"GTC"[..], 5).into(),
            (&b"AAA"[..]).into()
        ];

        let model = Hmm::from(&modules);
        let expected: Array2<f32> = read_npy("data/test/log_trans_f32.npy").unwrap();
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
        let initial = read_npy("data/test/log_init_f32.npy").unwrap();
        let transition = read_npy("data/test/log_trans_f32.npy").unwrap();
        let emission = read_npy("data/test/log_emit_f32.npy").unwrap();

        let model = Hmm { 
            states: Vec::new(),
            state_to_mod: Vec::new(),
            deletions: HashSet::new(), module_changes: HashSet::new(), unit_changes: HashSet::new(),
            initial, transition, emission 
        };
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
        use State as S;
        let states = vec![
            S::Start,
            S::Seq{c: b'T', id: 0}, S::Seq{c: b'C', id: 0}, S::Seq{c: b'T', id: 0},
            S::Motif{c: b'G', id: 1}, S::Motif{c: b'T', id: 0}, S::Motif{c: b'C', id: 1},
            S::Seq{c: b'A', id: 2}, S::Seq{c: b'A', id: 2}, S::Seq{c: b'A', id: 2},
            S::End,
            S::Ins, S::Ins, S::Ins, S::Ins,
            S::Ins, S::Ins, S::Ins, S::Ins
        ];
        let obtained = emission_probabilities(&states).map(|&x| x.ln());
        let expected: Array3<f32> = read_npy("data/test/log_emit_f32.npy").unwrap();
        assert_eq_ndarray3(expected.view(), obtained.view(), (1e-3, 2));
    }

    #[test]
    fn emissions_support_iupac_codes() {
        let modules: Vec<Module> = vec![
            (&b"TCTTGCTACG"[..]).into(),
            // (&b"GCA"[..], 5).into(),
            // (&b"GCR"[..], 5).into(),
            (&b"GCN"[..], 5).into(),
            (&b"TTCCCGGCTA"[..]).into()
        ];

        let model = Hmm::from(&modules);

        let seq =  b"TAGCTCTTGCTACGGCGGCAGCGGCAGCGGCAGCATTCCCGGCTATGT".to_vec();
        let qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec();

        let (_, annotation) = model.log_predict(&seq, &qual);

        let (new_annot, read) = model.realign(&annotation, &seq); 
        let rref = model.reconstruct_sequence(&new_annot);
        let mods = model.reconstruct_mod_ids(&new_annot);

        use std::str;
        println!("{}", str::from_utf8(&read).unwrap());
        println!("{}", str::from_utf8(&rref).unwrap());
        println!("{}", str::from_utf8(&mods).unwrap());
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

    // use ndarray::Dimension;
    // use ndarray::RemoveAxis;
    // maybe look at how outer_iter is implemented?
    // fn find_diff<D>(
    //     a1: ArrayView<f32, D>, a2: ArrayView<f32, D>, acc: (f32, i32)
    // ) -> Option<Vec<usize>>
    // where
    //     D: Dimension + RemoveAxis
    // {
    //     let tmp = a1.outer_iter();
    //     let shp = a1.shape();
    //     if shp.len() == 1 {
    //         for i in a1.outer_iter() {
    //         }
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

