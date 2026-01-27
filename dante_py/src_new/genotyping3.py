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
    max_rep = max(obs_counts) + 5

    model = Model(read_lengths, max_rep)
    likelihoods = model.evaluate(obs_counts, read_lengths, is_spanning, monoallelic_motif)

    _min_rep = get_min_rep(spanning_observed_counts)
    _max_rep = get_max_rep(spanning_observed_counts)
    likelihoods = transform_to_old_format(likelihoods, _min_rep, _max_rep, model.exp_idx, model.bkg_idx)
    predicted_tmp = predict(likelihoods)
    raw_confidence = get_confidence(likelihoods, predicted_tmp, _max_rep, monoallelic_motif)
    prediction = convert_to_sym(_max_rep, predicted_tmp, monoallelic_motif)

    return likelihoods, prediction, raw_confidence


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
    L_BCKG_MODEL_OPEN = 0.5

    # l_read_given_one_genotype and l_read_given_two_genotypes
    L_BCKG_OPEN = 0.01
    L_BCKG_CLOSED = 0.001  # L_BCKG_OPEN / 10

    # model_full
    P_DEL1, P_DEL2 = 0.0001, 0.0001
    P_INS = 0.0001

    def __init__(self, read_lengths: list[int], max_rep: int):
        """ Initialization of the Inference class + setup of all models and their probabilities. """
        self.read_dist: np.ndarray = np.bincount(read_lengths, minlength=100) / len(read_lengths)

        # max_spanning_rep
        # max_flanking_rep

        self.max_rep: int = max_rep
        self.exp_idx: int = max_rep + 1
        self.bkg_idx: int = max_rep + 2
        self.mprobs: list = Model.construct_mprobs(self.max_rep)    # P(G)
        self.models: list = Model.construct_models(self.max_rep)    # P(A|G)

    @staticmethod
    def construct_mprobs(max_rep: int) -> list[float]:
        """ Construct model probabilities (float - not summing to 1? because they are likelihoods?) """
        mprobs = []
        for i in range(max_rep + 1):
            mprobs.append(Model.L_OTHERS)
        mprobs.append(Model.L_EXP)
        mprobs.append(Model.L_BCKG_MODEL_OPEN)
        return mprobs

    @staticmethod
    def construct_models(max_rep: int) -> list[np.ndarray]:
        models = []
        for i in range(max_rep + 1):  # inclusive 0, 1, ..., n
            models.append(Model.model_full(max_rep, i))
        models.append(Model.model_full(max_rep, max_rep - 1))   # exp  this should be something else
        models.append(Model.model_bckg(max_rep))                # bkg
        return models

    @staticmethod
    def model_full(size: int, gt: int) -> np.ndarray:
        """Returns ndarray with length size"""
        def clip(value, minimal, maximal):
            return min(max(minimal, value), maximal)

        p_del = clip(Model.P_DEL1 + Model.P_DEL2 * gt, 0.0, 1.0)
        deletes = binom.pmf(np.arange(gt + 1), gt, p_del)  # this should be geometric distribution
        p_ins = Model.P_INS
        inserts = binom.pmf(np.arange(gt + 1), gt, p_ins)

        result = np.convolve(inserts, deletes[::-1])[:size]
        padding = np.zeros(size - len(result), dtype=float)
        return np.concatenate([result, padding])

    @staticmethod
    def model_bckg(size: int) -> np.ndarray:
        """Returns ndarray with length size"""
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
            bkground_likelihood = self.L_BCKG_CLOSED if sf else self.L_BCKG_OPEN
            bg_idx = self.bkg_idx
            bckgrnd_l = bkground_likelihood * self.l_read_given_genotype(oc, rl, sf, bg_idx)
            allele1_l = self.mprobs[g1_idx] * self.l_read_given_genotype(oc, rl, sf, g1_idx)
            allele2_l = self.mprobs[g2_idx] * self.l_read_given_genotype(oc, rl, sf, g2_idx)
            value = bckgrnd_l + allele1_l + allele2_l
            m_lh += np.log(value)
        return m_lh

    @functools.lru_cache()
    def l_read_given_genotype(self, oc: int, rl: int, is_spanning: bool, gt_idx: int) -> float:
        if is_spanning:
            return self.l_spanning_read_given_genotype(oc, rl, gt_idx)
        else:
            return self.l_flanking_read_given_genotype(oc, rl, gt_idx)

    def l_spanning_read_given_genotype(self, oc: int, rl: int, gt_idx: int) -> float:
        """ This wants to be eq. 6 in https://doi.org/10.1093/bioinformatics/bty791 """
        likelihood_rl: float = self.read_dist[rl]
        likelihood_model: float = self.models[gt_idx][oc]
        likelihood_cov: float = 1.0 / rl
        return likelihood_rl * likelihood_model * likelihood_cov

    def l_flanking_read_given_genotype(self, oc: int, rl: int, gt_idx: int) -> float:
        """ This wants to be eq. 6 in https://doi.org/10.1093/bioinformatics/bty791 for flanking reads"""
        likelihood_rl: float = self.read_dist[rl]

        partial_likelihood = 0.0
        number_of_options = 0
        for true_length in list(range(oc, self.max_rep)):
            likelihood_model: float = self.models[gt_idx][true_length]
            likelihood_cov: float = 1.0 / rl

            partial_likelihood += likelihood_model * likelihood_cov
            number_of_options += 1

        return likelihood_rl * partial_likelihood / float(number_of_options)
# ---


def transform_to_old_format(lhoods, min_rep, max_rep, exp_idx, bkg_idx):
    likelihoods = np.zeros((max_rep, max_rep + 1))
    rng = slice(min_rep, max_rep)
    likelihoods[rng, rng] = lhoods[rng, rng]
    likelihoods[0, 0] = lhoods[bkg_idx, bkg_idx]
    likelihoods[0, max_rep] = lhoods[exp_idx, exp_idx]
    likelihoods[rng, max_rep] = lhoods[rng, exp_idx]
    likelihoods[0, rng] = lhoods[rng, bkg_idx]
    ind_good = (likelihoods < 0.0) & (likelihoods > -1e10) & (likelihoods != np.nan)
    likelihoods[~ind_good] = -np.inf
    return likelihoods


def predict(lh_array: np.ndarray) -> tuple[int, int]:
    ind_good = lh_array != -np.inf
    if len(lh_array[ind_good]) == 0:
        return 0, 0
    best = sorted(np.unravel_index(np.argmax(lh_array), lh_array.shape))
    prediction = (int(best[0]), int(best[1]))

    return prediction


def convert_to_sym(max_rep, best: tuple[int, int], monoallelic: bool) -> tuple[int | str, int | str]:
    """ Convert numeric alleles to their symbolic representations. """
    def fn1(x):
        return 'E' if x == max_rep else 'B' if x == 0 else x

    if best[0] == 0 and best[1] == max_rep:
        best_sym = ('E', 'E')
    else:
        best_sym = tuple(map(fn1, best))

    # if mono-allelic return 'X' as second allele symbol
    if monoallelic:
        best_sym = (best_sym[0], 'X')

    return best_sym


def get_confidence(lh_array: np.ndarray, predicted: tuple[int, int], max_rep: int, monoallelic: bool) -> Confidences:
    """ Get confidence of a prediction. """
    lh_corr_array = lh_array - np.max(lh_array)
    lh_sum = np.sum(np.exp(lh_corr_array))

    confidence: float = np.exp(lh_corr_array[predicted[0], predicted[1]]) / lh_sum
    confidence1: float
    confidence2: float
    if predicted[0] == predicted[1]:  # same alleles - we compute the probability per allele
        confidence1 = np.sum(np.exp(lh_corr_array[predicted[0], :])) / lh_sum
        confidence2 = np.sum(np.exp(lh_corr_array[:, predicted[1]])) / lh_sum
    elif predicted[1] == lh_corr_array.shape[0]:  # expanded allele - expanded is only on one side of the array
        confidence1 = (
            np.sum(np.exp(lh_corr_array[predicted[0], :]))
            + np.sum(np.exp(lh_corr_array[:, predicted[0]]))
            - np.exp(lh_corr_array[predicted[0], predicted[0]])
        ) / lh_sum
        confidence2 = np.sum(np.exp(lh_corr_array[:, predicted[1]])) / lh_sum
    else:  # normal behavior - different alleles , no expanded, compute all likelihoods of the alleles
        confidence1 = (
            np.sum(np.exp(lh_corr_array[predicted[0], :]))
            + np.sum(np.exp(lh_corr_array[:, predicted[0]]))
            - np.exp(lh_corr_array[predicted[0], predicted[0]])
        ) / lh_sum
        confidence2 = (
            np.sum(np.exp(lh_corr_array[:, predicted[1]]))
            + np.sum(np.exp(lh_corr_array[predicted[1], :]))
            - np.exp(lh_corr_array[predicted[1], predicted[1]])
        ) / lh_sum

    confidence_back: float = np.exp(lh_corr_array[0, 0]) / lh_sum
    confidence_back_all: float = np.sum(np.exp(lh_corr_array[0, :])) / lh_sum
    confidence_exp: float = np.exp(lh_corr_array[0, max_rep]) / lh_sum
    confidence_exp_all: float = np.sum(np.exp(lh_corr_array[:, max_rep])) / lh_sum

    if monoallelic:
        confidence2 = np.nan

    result = (
        confidence, confidence1, confidence2,
        confidence_back, confidence_back_all, confidence_exp, confidence_exp_all
    )
    return result


MIN_REPETITIONS = 1
OVERHEAD = 3


def get_min_rep(spanning_obs_counts: list[int]) -> int:
    return max(MIN_REPETITIONS, min(spanning_obs_counts) - OVERHEAD)  # inclusive


def get_max_rep(spanning_obs_counts: list[int]) -> int:
    return max(spanning_obs_counts) + OVERHEAD + 1  # non-inclusive
