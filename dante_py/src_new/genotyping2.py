"""module for genotypization"""

import functools
import itertools
from typing import Iterator, TypeAlias

import numpy as np
from scipy.stats import binom  # type: ignore
np.set_printoptions(precision=4, suppress=True, linewidth=300, floatmode='fixed')

Confidences: TypeAlias = tuple[float, float, float, float, float, float, float]


def genotype(
    spanning_observed_counts: list[int], spanning_read_lengths: list[int],
    flanking_observed_counts: list[int], flanking_read_lengths: list[int],
    monoallelic_motif: bool, min_rep_count: int, min_flank_len: int, min_rep_len: int
) -> tuple:
    """This function provides an interface for genotypization step."""
    print()
    if len(spanning_observed_counts) == 0:
        return (None, ('B', 'B'), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

    model = Inference(spanning_observed_counts, spanning_read_lengths, flanking_observed_counts, flanking_read_lengths)

    likelihoods = model.evaluate(
        spanning_observed_counts, spanning_read_lengths,
        flanking_observed_counts, flanking_read_lengths,
        monoallelic_motif
    )

    predicted_tmp = predict(likelihoods)
    raw_confidence = get_confidence(likelihoods, predicted_tmp, model.max_rep, monoallelic_motif)
    prediction = convert_to_sym(model.max_rep, predicted_tmp, monoallelic_motif)

    return likelihoods, prediction, raw_confidence


# All other objects below this line are considered internal a should not be used
class Inference:
    """ Class for inference of alleles. """
    MIN_REPETITIONS = 1
    OVERHEAD = 3

    L_EXP = 1.01
    L_OTHERS = 1.0
    L_BCKG_OPEN = 0.01
    L_BCKG_CLOSED = 0.001  # L_BCKG_OPEN / 10

    L_BCKG_MODEL_OPEN = 0.5

    def __init__(
        self,
        spanning_obs_counts: list[int], spanning_read_lengths: list[int],
        flanking_obs_counts: list[int], flanking_read_lengths: list[int]
    ):
        """
        Initialization of the Inference class + setup of all models and their probabilities.
        """
        read_distribution = np.bincount(spanning_read_lengths + flanking_read_lengths, minlength=100)
        self.read_dist: np.ndarray = read_distribution / float(np.sum(read_distribution))  # make it sum to 1.0

        # TODO: Can I remove this and have 0, 1, ..., n, E, B?
        min_rep, max_rep, e_allele = self.get_boundaries(spanning_obs_counts, flanking_obs_counts)
        self.min_rep = min_rep
        self.max_with_e = e_allele

        self.max_rep = max_rep  # should be inclusive, but I think it often isn't
        self.exp_idx = max_rep + 1
        self.bkg_idx = max_rep + 2

        models, mprobs = Inference.construct_models(self.min_rep, self.max_rep, self.max_with_e)
        self.models = models
        self.mprobs = mprobs

    @classmethod
    def get_boundaries(cls, spanning_obs_counts: list[int], flanking_obs_counts: list[int]) -> tuple[int, int, int]:
        """Returns boundaries of the implicit matrix. Is it possible to get rid of this?"""
        max_rep = max(spanning_obs_counts) + cls.OVERHEAD + 1  # non-inclusive
        min_rep = max(cls.MIN_REPETITIONS, min(spanning_obs_counts) - cls.OVERHEAD)  # inclusive

        # expanded allele
        max_with_e = max_rep + 1
        if len(flanking_obs_counts) > 0:
            max_with_e = max(max_rep, max(flanking_obs_counts) + 1) + 1
        return (min_rep, max_rep, max_with_e)

    @staticmethod
    def construct_models(min_rep, max_rep, max_with_e) -> tuple[list, list]:
        """
        Construct models (np.ndarray representing probability distribution of getting read given haplotype) and
        model probabilities (float - not summing to 1? because they are likelihoods?)
        """
        # get new models
        new_models = []
        for i in range(max_rep + 1):  # inclusive 0, 1, ..., n
            new_models.append(model_full(max_with_e, i))
        new_models.append(model_full(max_with_e, max_with_e - 1))     # exp
        new_models.append(model_bckg(max_with_e, min_rep))            # bkg

        # get new mprobs
        new_mprobs = []
        for i in range(max_rep + 1):
            new_mprobs.append(Inference.L_OTHERS)
        new_mprobs.append(Inference.L_EXP)
        new_mprobs.append(Inference.L_BCKG_MODEL_OPEN)

        return new_models, new_mprobs

    def evaluate(
        self,
        spanning_obs_counts: list[int], spanning_read_lengths: list[int],
        flanking_obs_counts: list[int], flanking_read_lengths: list[int],
        is_monoallelic: bool
    ) -> np.ndarray:
        obs_counts = np.array(spanning_obs_counts + flanking_obs_counts)
        read_lengths = np.array(spanning_read_lengths + flanking_read_lengths)
        is_spanning = np.array([True] * len(spanning_obs_counts) + [False] * len(flanking_obs_counts))

        if is_monoallelic:
            return self.evaluate_monoallelic_motif(obs_counts, read_lengths, is_spanning)
        else:
            return self.evaluate_biallelic_motif(obs_counts, read_lengths, is_spanning)
    # ---

    def evaluate_monoallelic_motif(self, obs_counts, read_lengths, is_spanning):
        n = self.max_rep + 3  # 0, 1, ..., n, E, B
        lhoods = np.full((n, n), -np.inf)  # because they are loglikelihoods
        for idx in range(n):
            m_lh = 0.0
            for obs, rl, closed in zip(obs_counts, read_lengths, is_spanning):
                m_lh += np.log(self.l_read_given_one_genotype(obs, rl, closed, idx))
            lhoods[idx, idx] = m_lh

        # transform to old format
        likelihoods = np.zeros((self.max_rep, self.max_rep + 1))
        for idx in range(self.min_rep, self.max_rep):
            likelihoods[idx, idx] = lhoods[idx, idx]
        likelihoods[0, 0] = lhoods[self.bkg_idx, self.bkg_idx]
        likelihoods[0, self.max_rep] = lhoods[self.exp_idx, self.exp_idx]
        ind_good = (likelihoods < 0.0) & (likelihoods > -1e10) & (likelihoods != np.nan)
        likelihoods[~ind_good] = -np.inf
        return likelihoods

    def evaluate_biallelic_motif(self, obs_counts, read_lengths, is_spanning):
        """
        This description is wrong, but slightly useful.

        Evaluates all the models in explicit matrix.
        Does all the inference,
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
        def generate_biallelic_indices(min_rep: int, max_rep: int) -> Iterator:
            # B = 0
            # E = max_rep
            for model_index1 in range(min_rep, max_rep):
                for model_index2 in range(model_index1, max_rep):
                    yield (model_index1, model_index2)  # , (model_index1, model_index2)
                yield (model_index1, max_rep)  # , (model_index1, 'E')
                yield (0, model_index1)  # , ('B', model_index1)

            yield (0, 0)  # , ('B', 'B')
            yield (max_rep, max_rep)  # , ('E', 'E')

        likelihoods = np.zeros((self.max_rep, self.max_rep + 1))
        models = list(generate_biallelic_indices(self.min_rep, self.max_rep))
        # print(models)
        for m1, m2 in models:
            # P(OC, RL, SF | G1, G2), where OC=obs_counts, RL=read_lengths, SF=is_spanning, G1=m1, G2=m2
            model_lh = 0.0
            for obs, rl, closed in zip(obs_counts, read_lengths, is_spanning):
                # lh = self.likelihood_read(obs, rl, m1, m2, closed=closed)
                m1_new = self.bkg_idx if m1 == 0 else self.exp_idx if m1 == self.max_rep else m1
                m2_new = self.bkg_idx if m2 == 0 else self.exp_idx if m2 == self.max_rep else m2
                lh = self.l_read_given_two_genotypes(obs, rl, closed, m1_new, m2_new)
                model_lh += np.log(lh)

            if m1 == self.max_rep and m2 == self.max_rep:  # legacy, for backwards compatibility, change in future
                m1 = 0
            likelihoods[m1, m2] = model_lh

        ind_good = (likelihoods < 0.0) & (likelihoods > -1e10) & (likelihoods != np.nan)
        likelihoods[~ind_good] = -np.inf
        return likelihoods

    # ---
    @functools.lru_cache()
    def l_read_given_two_genotypes(
        self, oc: int, rl: int, sf: bool, g1_idx: int, g2_idx: int
    ) -> float:
        """ P(OC[i], RL[i], SF[i] | G1, G2) * P(G1, G2) which is definitelly incorrect"""
        bkground_likelihood = self.L_BCKG_CLOSED if sf else self.L_BCKG_OPEN
        bg_idx = self.bkg_idx
        bckgrnd_l = bkground_likelihood * self.l_read_given_genotype(oc, rl, sf, bg_idx)
        allele1_l = self.mprobs[g1_idx] * self.l_read_given_genotype(oc, rl, sf, g1_idx)
        allele2_l = self.mprobs[g2_idx] * self.l_read_given_genotype(oc, rl, sf, g2_idx)

        # TODO: tuto podla mna nemoze byt len tak +, chyba tam korelacia modelov, ale v ramci zjednodusenia asi ok
        return allele1_l + allele2_l + bckgrnd_l
        # return bckgrnd_l + allele1_l + allele2_l  # damn you float arithmetic

    @functools.lru_cache()
    def l_read_given_one_genotype(
        self, oc: int, rl: int, sf: bool, g1_idx: int
    ) -> float:
        """ P(OC[i], RL[i], SF[i] | G1) * P(G1) """
        # print(g1_idx, sep=" ", end="")
        bkground_likelihood = self.L_BCKG_CLOSED if sf else self.L_BCKG_OPEN
        bg_idx = self.bkg_idx
        bckgrnd_l = bkground_likelihood * self.l_read_given_genotype(oc, rl, sf, bg_idx)
        allele1_l = self.mprobs[g1_idx] * self.l_read_given_genotype(oc, rl, sf, g1_idx)

        # TODO: tuto podla mna nemoze byt len tak +, chyba tam korelacia modelov, ale v ramci zjednodusenia asi ok
        return bckgrnd_l + allele1_l

    def l_read_given_genotype(self, oc: int, rl: int, is_spanning: bool, gt_idx: int) -> float:
        if is_spanning:
            return self.l_spanning_read_given_genotype(oc, rl, gt_idx)
        else:
            return self.l_flanking_read_given_genotype(oc, rl, gt_idx)

    def l_spanning_read_given_genotype(self, oc, rl, gt_idx) -> float:
        """ This wants to be eq. 6 in https://doi.org/10.1093/bioinformatics/bty791 """

        def lc(oc, rl):  # basically returns 1/rl with some inherited confusion numbers
            return 1.0 / float(max(0, +(rl - 5 + oc)) - max(0, -(rl - 7 - oc)))  # This never made sense anyway

        likelihood_rl: float = self.read_dist[rl]
        likelihood_model: float = self.models[gt_idx][oc]
        likelihood_cov: float = lc(oc, rl)
        return likelihood_rl * likelihood_model * likelihood_cov

    def l_flanking_read_given_genotype(self, oc, rl, gt_idx) -> float:
        """ This wants to be eq. 6 in https://doi.org/10.1093/bioinformatics/bty791 for flanking reads"""

        def lc(oc, rl):  # basically returns 1/rl with some inherited confusion numbers
            return 1.0 / float(max(0, +(rl - 5 + oc)) - max(0, -(rl - 7 - oc)))  # This never made sense anyway

        likelihood_rl: float = self.read_dist[rl]

        partial_likelihood = 0.0
        number_of_options = 0
        tmp = list(itertools.chain(range(oc, self.max_rep), [self.max_with_e - 1]))
        for true_length in tmp:
            likelihood_model: float = self.models[gt_idx][true_length]  # I think this is supposed to be other way
            likelihood_cov: float = lc(true_length, rl)

            partial_likelihood += likelihood_model * likelihood_cov
            number_of_options += 1

        return likelihood_rl * partial_likelihood / float(number_of_options)
    # ---


def predict(lh_array: np.ndarray) -> tuple[int, int]:
    ind_good = lh_array != -np.inf
    if len(lh_array[ind_good]) == 0:
        return 0, 0
    best = sorted(np.unravel_index(np.argmax(lh_array), lh_array.shape))
    prediction = (int(best[0]), int(best[1]))

    return prediction


def convert_to_sym(max_rep, best: tuple[int, int], monoallelic: bool) -> tuple[int | str, int | str]:
    """
    Convert numeric alleles to their symbolic representations.
    :param best: (int, int) - numeric representation of alleles
    :param monoallelic: bool - if this is monoallelic version
    :return: (int|str, int|str) - symbolic representation of alleles
    """
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
    """
    Get confidence of a prediction.
    :param lh_array: 2D-ndarray - log likelihoods of the prediction
    :param predicted: tuple(int, int) - predicted alleles
    :param monoallelic: bool - do we have a mono-allelic motif (i.e. chrX/chrY and male sample?)
    :return: tuple[float, float, float | str, float, float, float, float] - prediction confidence of
    all, first, and second allele(s), background and expanded states
    """
    # get confidence
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
        # confidence2 = '---'  # TODO: fix this
        confidence2 = np.nan

    return (
        confidence, confidence1, confidence2,
        confidence_back, confidence_back_all,
        confidence_exp, confidence_exp_all
    )


# ---
def model_full(size: int, gt: int) -> np.ndarray:
    """
    Create binomial model for both deletes and inserts of STRs
    :param rng: int - max_range of distribution
    :param n: int - target allele number
    :return: ndarray - combined distribution
    """
    def clip(value, minimal, maximal):
        return min(max(minimal, value), maximal)

    def combine_distribs(deletes: np.ndarray, inserts: np.ndarray) -> np.ndarray:
        # how much to fill?
        to_fill = sum(deletes == 0.0) + 1
        while to_fill < len(inserts) and inserts[to_fill] > 0.0001:
            to_fill += 1

        # create the end array
        end_distr = np.zeros_like(deletes, dtype=float)

        # fill it!
        for i, a in enumerate(inserts[:to_fill]):
            end_distr[i:] += (deletes * a)[:len(deletes) - i]

        return end_distr

    # DEFAULT_MODEL_PARAMS = (0.001, 0.000105087, 0.0210812, 0.001)
    # p1, p2, p3, q = DEFAULT_MODEL_PARAMS
    # inserts = q
    # deletes = p1 + p2 * n
    DEFAULT_MODEL_PARAMS = (0.0001, 0.0001, 0.0, 0.0001)

    p1, p2, _, q = DEFAULT_MODEL_PARAMS
    deletes: np.ndarray = binom.pmf(np.arange(size), gt, clip(1 - (p1 + p2 * gt), 0.0, 1.0))
    # print(deletes)
    inserts: np.ndarray = binom.pmf(np.arange(size), gt, q)
    # print(inserts)
    result = combine_distribs(deletes, inserts)
    # print(result)
    return result


def model_bckg(size: int, zeros: int) -> np.ndarray:
    """Returns ndarray with length size"""
    result = np.concatenate([
        np.zeros(zeros, dtype=float),
        np.ones(size - zeros, dtype=float) / float(size - zeros)
    ])
    return result
