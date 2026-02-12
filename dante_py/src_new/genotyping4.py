"""module for genotypization"""

from scipy.stats import binom  # type: ignore
from typing import TypeAlias
import functools
import numpy as np
from typing import Iterable

Confidences: TypeAlias = tuple[float, float, float, float, float, float, float]


def genotype(
    spanning_observed_counts: list[int], spanning_read_lengths: list[int],
    flanking_observed_counts: list[int], flanking_read_lengths: list[int],
    monoallelic_motif: bool, _min_rep_count: int, _min_flank_len: int, _min_rep_len: int
) -> tuple:
    """This function provides an interface for genotypization step."""
    print()
    if len(spanning_observed_counts) == 0:
        return (None, ('B', 'B'), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

    obs_counts = spanning_observed_counts + flanking_observed_counts
    read_lengths = spanning_read_lengths + flanking_read_lengths
    is_spanning = [True] * len(spanning_observed_counts) + [False] * len(flanking_observed_counts)
    max_spanning_reps = max(spanning_observed_counts)
    max_overall_reps = max(obs_counts)

    model = Model(read_lengths, max_spanning_reps, max_overall_reps)
    likelihoods = model.evaluate(obs_counts, read_lengths, is_spanning, monoallelic_motif)
    pred_sym = model.predict_sym(likelihoods, monoallelic_motif)
    confidences = model.get_conf(likelihoods, monoallelic_motif)

    _max_rep = model.max_rep + 1
    _min_rep = get_min_rep(spanning_observed_counts)
    likelihoods = transform_to_old_format(likelihoods, _min_rep, _max_rep, model.exp_idx, model.bkg_idx)

    return likelihoods, pred_sym, confidences


# All other objects below this line are considered internal a should not be used
class Model:
    """
    Class for inference of alleles.
    This description is wrong, but slightly useful.

    computes for which 2 combination of alleles are these annotations and parameters the best.
    argmax_{G1, G2} P(G1, G2 | AL, COV, RL)
        ~ P(AL, COV, RL | G1, G2) * P(G1, G2)
        = prod_{read_i} P(al_i, cov_i, rl_i | G1, G2) * P(G1, G2)
        = independent G1 G2
        = prod_{read_i} P(al_i, cov_i, rl_i | G1) * P(al_i, cov_i, rl_i | G2) * P(G1) * P(G2)
        {here G1, G2 is from possible alleles, background, and expanded, priors are from params}

     P(al_i, cov_i, rl_i | G1) - 2 options:
         1. closed evidence (al_i = X), we know X;
         2. open evidence (al_i >= X), cl_i == True if i is closed

     1.: P(al_i, cov_i, rl_i, cl_i | G1)
        = P(rl_i is from read distrib.) * p(allele is == al_i | G1) * P(read generated closed evidence | rl_i, al_i)
     2.: P(al_i, cov_i, rl_i, cl_i | G1)
        = P(rl_i is from read distrib.) * P(allele is >= al_i | G1) * P(read generated open evidence | rl_i, al_i)
    """
    # construct_mprobs
    L_OTHERS = 1.0
    L_EXP = 1.01
    L_BKG = 0.01

    # model_full
    P_DEL1, P_DEL2 = 0.0001, 0.0001
    P_INS = 0.0001

    def __init__(self, read_lengths: list[int], max_spanning_rep: int, max_flanking_rep: int):
        """ Initialization of the Inference class + setup of all models and their probabilities. """
        self.read_dist: np.ndarray = np.bincount(read_lengths, minlength=100) / len(read_lengths)
        self.max_rep: int = max_spanning_rep
        self.max_frep: int = max_flanking_rep

        self.exp_idx: int = self.max_rep + 1
        self.bkg_idx: int = self.max_rep + 2
        self.mprobs: list = Model.construct_mprobs(self.max_rep)                    # P(G)
        self.models: list = Model.construct_models(self.max_rep, self.max_frep)     # P(A|G)

    @staticmethod
    def construct_mprobs(max_rep: int) -> list[float]:
        """ Construct model probabilities (float - not summing to 1? because they are likelihoods?) """
        mprobs = []
        for i in range(max_rep + 1):
            mprobs.append(Model.L_OTHERS)
        mprobs.append(Model.L_EXP)
        mprobs.append(Model.L_BKG)
        # mprobs should sum to 1.0
        # should be either uniform or from https://gnomad.broadinstitute.org/short-tandem-repeat/ATXN1
        return mprobs

    @staticmethod
    def construct_models(max_rep: int, max_flanking_rep: int) -> list[np.ndarray]:
        models = []
        for i in range(max_rep + 1):  # inclusive 0, 1, ..., n
            models.append(Model.model_full(max_flanking_rep + 1, i))
        models.append(Model.model_expn(max_flanking_rep + 1, max_rep + 1))  # exp
        models.append(Model.model_bckg(max_flanking_rep + 1))               # bkg
        return models

    @staticmethod
    def model_full(size: int, gt: int) -> np.ndarray:
        """Returns ndarray with length size"""
        def clip(value, minimal, maximal):
            return min(max(minimal, value), maximal)

        p_del = clip(Model.P_DEL1 + Model.P_DEL2 * gt, 0.0, 1.0)
        deletes = binom.pmf(np.arange(gt + 1), gt, p_del)  # this should be geometric distribution
        p_ins = Model.P_INS
        inserts = binom.pmf(np.arange(gt + 1), gt, p_ins)  # this should be geometric distribution

        result = np.convolve(inserts, deletes[::-1])[:size]
        padding = np.zeros(size - len(result), dtype=float)
        return np.concatenate([result, padding])

    @staticmethod
    def model_expn(size: int, gt_min: int) -> np.ndarray:
        """Returns ndarray with length size"""
        # expansion should be sum i from (n+1) to (inf) model_full(n, i)
        result = np.zeros(size, dtype=float)
        for i in range(gt_min, size + 1):
            result += Model.model_full(size, i)
        result /= (size - gt_min + 1)
        return result

    @staticmethod
    def model_bckg(size: int) -> np.ndarray:
        """Returns ndarray with length size"""
        # representing microsatellite instability?
        return np.ones(size, dtype=float) / float(size)
    # ---

    def evaluate(
        self, observed: list[int], rlengths: list[int], spanning: list[bool], is_monoallelic: bool
    ) -> np.ndarray:
        """Returns a matrix of loglikelihoods for each considered option"""
        n = self.max_rep + 3  # 0, 1, ..., n, E, B
        llmatrix = np.full((n, n), -np.inf)
        genotypes: Iterable
        if is_monoallelic:
            genotypes = range(n)
            for gt_idx in genotypes:
                llmatrix[gt_idx, gt_idx] = self.loglikelihood_of_D_given_G(observed, rlengths, spanning, gt_idx, gt_idx)
        else:
            genotypes = ((g1_idx, g2_idx) for g1_idx in range(n) for g2_idx in range(g1_idx, n))
            for (g1_idx, g2_idx) in genotypes:
                llmatrix[g1_idx, g2_idx] = self.loglikelihood_of_D_given_G(observed, rlengths, spanning, g1_idx, g2_idx)
        return llmatrix

    def loglikelihood_of_D_given_G(
        self, obs_counts: list[int], read_lengths: list[int], is_spanning: list[bool], g1_idx: int, g2_idx: int
    ) -> float:
        """ This wants to be eq. 6 in https://doi.org/10.1093/bioinformatics/bty791 """
        """ P(OC, RL, SF | G1, G2) """
        m_lh = 0.0
        for oc, rl, sf in zip(obs_counts, read_lengths, is_spanning):
            bckgrnd_l = self.l_read_given_genotype(oc, rl, sf, self.bkg_idx)
            allele1_l = self.l_read_given_genotype(oc, rl, sf, g1_idx)
            allele2_l = self.l_read_given_genotype(oc, rl, sf, g2_idx)
            m_lh += np.log(bckgrnd_l + allele1_l + allele2_l)
        return m_lh

    @functools.lru_cache()
    def l_read_given_genotype(self, oc: int, rl: int, is_spanning: bool, gt_idx: int) -> float:
        """ This wants to be eq. 6 in https://doi.org/10.1093/bioinformatics/bty791 """
        """ P(oc, rl, sf | G) """
        lh_cover: float = 1.0 / rl              # P(b_i | a_i, r_i) # incorrect, but does something
        lh_r_len: float = self.read_dist[rl]    # P(r_i)            # correct, but does nothing
        lh_model: float                         # P(a_i | g_i)
        lh_mprob: float = self.mprobs[gt_idx]   # P(g_i)

        if is_spanning:
            lh_model = self.models[gt_idx][oc]
        else:
            lh_model = sum(self.models[gt_idx][oc:]) / len(self.models[gt_idx][oc:])

        # print(f"{lh_cover:.4f} {lh_model:.4f} {lh_mprob:.4f} {lh_r_len:.4f}")
        return lh_cover * lh_r_len * lh_model * lh_mprob
    # ---

    @staticmethod
    def predict(llmatrix: np.ndarray) -> tuple[int, int]:
        best = sorted(np.unravel_index(np.argmax(llmatrix), llmatrix.shape))
        prediction = (int(best[0]), int(best[1]))
        return prediction

    def predict_sym(self, llmatrix: np.ndarray, is_monoallelic: bool) -> tuple[int | str, int | str]:
        prediction = self.predict(llmatrix)
        sym_prediction = (
            'E' if prediction[0] == self.exp_idx else 'B' if prediction[0] == self.bkg_idx else prediction[0],
            'E' if prediction[1] == self.exp_idx else 'B' if prediction[1] == self.bkg_idx else prediction[1],
        )
        if is_monoallelic:
            sym_prediction = (sym_prediction[0], 'X')
        return sym_prediction

    def get_conf(self, llmatrix: np.ndarray, is_monoallelic: bool) -> tuple[float, ...]:
        """ Returns seven floats representing different confidences """
        # This trick is needed, because exp of large negative is zero in floats.
        llmatrix = llmatrix - np.max(llmatrix)
        # print(llmatrix / np.log(10.0))  # log-ratios in base 10
        # llmatrix = llmatrix / 100       # temperature scaling
        prob = np.exp(llmatrix)/np.sum(np.exp(llmatrix))  # softmax
        # print(prob)

        pred = self.predict(llmatrix)
        conf_pred: float = prob[pred[0], pred[1]]
        conf_al_1: float = np.sum(prob[pred[0], :]) + np.sum(prob[:, pred[0]]) - prob[pred[0], pred[0]]
        conf_al_2: float = np.sum(prob[pred[1], :]) + np.sum(prob[:, pred[1]]) - prob[pred[1], pred[1]]
        if is_monoallelic:
            conf_al_2 = float("nan")

        bkg = self.bkg_idx
        conf_bckg: float = prob[bkg, bkg]
        conf_bg_t: float = np.sum(prob[bkg, :]) + np.sum(prob[:, bkg]) - prob[bkg, bkg]

        exp = self.exp_idx
        conf_expn: float = prob[exp, exp]
        conf_ex_t: float = np.sum(prob[exp, :]) + np.sum(prob[:, exp]) - prob[exp, exp]

        return (conf_pred, conf_al_1, conf_al_2, conf_bckg, conf_bg_t, conf_expn, conf_ex_t)
# ---


MIN_REPETITIONS = 1
OVERHEAD = 3


def get_min_rep(spanning_obs_counts: list[int]) -> int:
    return max(MIN_REPETITIONS, min(spanning_obs_counts) - OVERHEAD)  # inclusive


# TODO: split this into somethings integratable to class and conversion to old
def transform_to_old_format(lhoods, min_rep, max_rep, exp_idx, bkg_idx):
    print(lhoods.shape)
    likelihoods = np.zeros((max_rep, max_rep + 1))
    rng = slice(min_rep, max_rep)
    # rng = slice(min_rep, exp_idx)
    likelihoods[rng, rng] = lhoods[rng, rng]
    likelihoods[0, 0] = lhoods[bkg_idx, bkg_idx]
    likelihoods[0, max_rep] = lhoods[exp_idx, exp_idx]
    likelihoods[rng, max_rep] = lhoods[rng, exp_idx]
    likelihoods[0, rng] = lhoods[rng, bkg_idx]
    ind_good = (likelihoods < 0.0) & (likelihoods > -1e10) & (likelihoods != np.nan)
    likelihoods[~ind_good] = -np.inf
    print(likelihoods.shape)
    return likelihoods
