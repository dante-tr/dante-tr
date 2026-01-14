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

    if len(spanning_observed_counts) == 0 and len(flanking_observed_counts) == 0:
        return (None, ('B', 'B'), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

    model = Inference(
        spanning_observed_counts, spanning_read_lengths,
        flanking_observed_counts, flanking_read_lengths,
        min_rep_count, min_flank_len, min_flank_len, min_rep_len
    )

    likelihoods = model.evaluate(
        spanning_observed_counts, spanning_read_lengths,
        flanking_observed_counts, flanking_read_lengths,
        monoallelic_motif
    )

    predicted_tmp = model.real_predict(likelihoods)
    # adjust for no spanning reads (should output Background)
    if len(spanning_observed_counts) == 0:
        predicted_tmp = (0, 0)
    prediction = convert_to_sym(model.max_rep, predicted_tmp, monoallelic_motif)

    raw_confidence = get_confidence(likelihoods, predicted_tmp, model.max_rep, monoallelic_motif)
    return likelihoods, prediction, raw_confidence


# All other objects below this line are considered internal a should not be used
class Inference:
    """ Class for inference of alleles. """
    MIN_REPETITIONS = 1
    OVERHEAD = 3

    OPEN_TO_CLOSED = 10.0
    L_OTHERS = 1.0
    L_BCKG_OPEN = 0.01
    L_EXP = 1.01
    L_BCKG_MODEL_OPEN = 0.5

    def __init__(
        self,
        observed_annots: list[int], spanning_read_lengths: list[int],
        observed_fa: list[int], flanking_read_lengths: list[int],
        min_rep_count: int, min_lflank_len: int, min_rflank_len: int, min_rep_len: int
    ):
        """
        Initialization of the Inference class + setup of all models and their probabilities.
        """
        # assign variables
        self.str_rep: int = min_rep_count
        self.minl_primer1: int = min_lflank_len
        self.minl_primer2: int = min_rflank_len
        self.minl_str: int = min_rep_len

        read_distribution = np.bincount(spanning_read_lengths + flanking_read_lengths, minlength=100)
        self.read_dist: np.ndarray = read_distribution / float(np.sum(read_distribution))  # make it sum to 1.0

        # TODO: Can I remove this and have 0, 1, ..., n, E, B?
        boundaries = self.get_boundaries(observed_annots, observed_fa)
        self.min_rep = boundaries[0]
        self.max_rep = boundaries[1]
        self.e_allele = boundaries[2]
        self.max_with_e = self.e_allele + 1  # non-inclusive

        tmp = self.construct_models()
        self.models = tmp[0]
        self.model_probabilities = tmp[1]

    @classmethod
    def get_boundaries(cls, observed_annots, observed_fa) -> tuple[int, int, int]:
        """Returns boundaries of the implicit matrix. Is it possible to get rid of this?"""
        if len(observed_annots) > 0:
            max_rep = max(observed_annots) + cls.OVERHEAD + 1  # non-inclusive
            min_rep = max(cls.MIN_REPETITIONS, min(observed_annots) - cls.OVERHEAD)  # inclusive
        else:
            max_rep = max(observed_fa) + cls.OVERHEAD  # non-inclusive
            min_rep = max(cls.MIN_REPETITIONS, max(observed_fa) - cls.OVERHEAD)  # inclusive

        # expanded allele
        if len(observed_fa) > 0:
            e_allele = max(max_rep, max(observed_fa) + 1)
        else:
            e_allele = max_rep
        return (min_rep, max_rep, e_allele)

    def construct_models(self) -> tuple[dict, dict]:
        """
        Construct models (np.ndarray representing probability distribution of getting read given haplotype) and
        model probabilities (float - not summing to 1? because they are likelihoods?)
        """
        # get models
        background_model = model_bckg(self.min_rep, self.max_with_e)
        expanded_model = model_full(self.max_with_e, self.max_with_e - 1)
        allele_models = {
            i: model_full(self.max_with_e, i) for i in range(self.min_rep, self.max_rep)
        }
        models = {
            'E': expanded_model,
            'B': background_model
        }
        models.update(allele_models)  # type: ignore

        # get model likelihoods
        allele_model_probabilities = {
            i: self.L_OTHERS for i in range(self.min_rep, self.max_rep)
        }
        model_probabilities = {
            'E': self.L_EXP,
            'B': self.L_BCKG_MODEL_OPEN
        }
        model_probabilities.update(allele_model_probabilities)  # type: ignore
        return models, model_probabilities

    # ---
    def evaluate(
        self,
        spanning_observed_counts: list[int], spanning_read_lengths: list[int],
        flanking_observed_counts: list[int], flanking_read_lengths: list[int],
        monoallelic_motif: bool
    ) -> np.ndarray:
        """Evaluates all the models in implicit matrix"""
        flag_spanning_flanking = np.concatenate([
            np.ones_like(spanning_observed_counts, dtype=bool),
            np.zeros_like(flanking_observed_counts, dtype=bool)
        ]).astype(bool)
        lh_dict = self.real_infer(
            monoallelic_motif,
            np.array(spanning_observed_counts + flanking_observed_counts),
            np.array(spanning_read_lengths + flanking_read_lengths),
            flag_spanning_flanking
        )
        likelihoods = self.convert_to_ndarray(lh_dict)
        return likelihoods

    def real_infer(
        self, monoallelic: bool, observed_arr: np.ndarray, rl_arr: np.ndarray, closed_arr: np.ndarray
    ) -> dict[tuple[int | str, int | str], float]:
        """
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
            = P(rl_i from read distrib.) * p(allele is al_i | G1) * P(read generated closed evidence | rl_i, al_i)
         2.: P(rl_i is from r.distr.) * P(allele is >= al_i | G1) * P(read generated open evidence | rl_i, al_i)

        :param annotations: list(Annotation) - closed annotated reads (both primers set)
        :param filt_annotations: list(Annotation) - open annotated reads (only one primer set)
        :param index_rep: int - index of a repetition
        :param verbose: bool - print more stuff?
        :param monoallelic: bool - do we have a mono-allelic motif (i.e. chrX/chrY and male sample?)
        :return: dict(tuple(int, int):float) - directory of model indices to their likelihood
        """
        # go through every model and evaluate:
        evaluated_models = {}
        if monoallelic:
            models = list(generate_models_one_allele(self.min_rep, self.max_rep))
        else:
            models = list(generate_models(self.min_rep, self.max_rep, multiple_bckgs=True))

        # print(models)
        for m1, m2 in models:
            evaluated_models[(m1, m2)] = 0.0
            # go through every read
            for obs, rl, closed in zip(observed_arr, rl_arr, closed_arr):
                lh = self.likelihood_read(obs, rl, m1, None if m2 == 'X' else m2, closed=closed)
                # TODO weighted sum according to the closeness/openness of reads?
                evaluated_models[(m1, m2)] += np.log(lh)

        return evaluated_models

    def convert_to_ndarray(self, lh_dict: dict[tuple[int | str, int | str], float]) -> np.ndarray:
        """Converts dictionary into matrix and selects the best prediction"""
        # convert to a numpy array:
        lh_array = np.zeros((self.max_rep, self.max_rep + 1))
        for (k1, k2), v in lh_dict.items():
            if k2 == 'X':  # if we have mono-allelic
                k2 = k1
            # B is the smallest, E is the largest!
            if k2 == 'B' or k1 == 'E' or (isinstance(k1, int) and isinstance(k2, int) and k2 < k1):
                k1, k2 = k2, k1
            if k1 == 'B':
                k1 = 0
            if k2 == 'B':
                k2 = 0
            if k1 == 'E':  # only if k2 is 'E' too.
                k1 = 0
            if k2 == 'E':
                k2 = self.max_rep
            lh_array[k1, k2] = v

        # get minimal and maximal likelihood
        ind_good = (lh_array < 0.0) & (lh_array > -1e10) & (lh_array != np.nan)
        if len(lh_array[ind_good]) == 0:
            return lh_array
        lh_array[~ind_good] = -np.inf

        # output best option
        return lh_array

    # TODO: remove this
    def likelihood_coverage(self, true_length, rl, _closed):
        """
        Likelihood of generating a read with this length and this allele.
        :param true_length: int - true number of repetitions of an STR
        :param rl: int - read length
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return: float - likelihood of a read being generated with this attributes
        """
        whole_inside_str = max(0, true_length * self.str_rep + self.minl_primer1 + self.minl_primer2 - rl + 1)
        # closed_overlapping = max(0, rl - self.minl_primer1 - self.minl_primer2 - true_length * self.str_rep + 1)
        open_overlapping = max(0, rl + true_length * self.str_rep - 2 * self.minl_str + 1)

        assert open_overlapping > whole_inside_str, \
            f"{open_overlapping} open {whole_inside_str} whole inside {true_length} {rl} {self.minl_str}"

        return 1.0 / float(open_overlapping - whole_inside_str)

    def likelihood_read_allele(self, model, observed, rl, closed=True) -> float:
        """
        Likelihood of generation of read with observed allele count and rl.
        :param model: ndarray - model for the allele
        :param observed: int - observed allele count
        :param rl: int - read length
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return:
        """
        if closed:
            likelihood_rl: float = self.read_dist[rl]
            likelihood_model = model[observed]
            likelihood_coverage = self.likelihood_coverage(observed, rl, True)
            return likelihood_rl * likelihood_model * likelihood_coverage

        number_of_options = 0
        partial_likelihood = 0
        for true_length in itertools.chain(range(observed, self.max_rep), [self.max_with_e - 1]):
            likelihood_model = model[true_length]
            likelihood_coverage = self.likelihood_coverage(true_length, rl, False)
            partial_likelihood += likelihood_model * likelihood_coverage
            number_of_options += 1

        likelihood_rl = self.read_dist[rl]
        return likelihood_rl * partial_likelihood / float(number_of_options)

    @functools.lru_cache()
    def likelihood_read(
        self, observed: int, rl: int, model_index1: int, model_index2: int | None = None, closed: bool = True
    ) -> float:
        """
        Compute likelihood of generation of a read from either of those models.
        :param observed: int - observed allele count
        :param rl: int - read length
        :param model_index1: char/int - model index for left allele
        :param model_index2: char/int - model index for right allele or None if mono-allelic
        :param closed: bool - if the read is closed - i.e. both primers are there
        :return: float - likelihood of this read generation
        """
        # TODO: tuto podla mna nemoze byt len tak +, chyba tam korelacia modelov, ale v ramci zjednodusenia asi ok
        m1 = model_index1
        m2 = model_index2
        allele1_likelihood = (
            self.model_probabilities[m1] * self.likelihood_read_allele(self.models[m1], observed, rl, closed)
        )
        allele2_likelihood = 0.0 if model_index2 is None else (
            self.model_probabilities[m2] * self.likelihood_read_allele(self.models[m2], observed, rl, closed)
        )

        if closed:
            p_bckg = self.L_BCKG_OPEN / self.OPEN_TO_CLOSED
        else:
            p_bckg = self.L_BCKG_OPEN
        bckgrnd_likelihood = p_bckg * self.likelihood_read_allele(self.models['B'], observed, rl, closed)

        assert not np.isnan(allele2_likelihood)
        assert not np.isnan(allele1_likelihood)
        assert not np.isnan(bckgrnd_likelihood)

        # "tuto" refers to the next line
        return allele1_likelihood + allele2_likelihood + bckgrnd_likelihood

    def real_predict(self, lh_array: np.ndarray) -> tuple[int, int]:
        # get minimal and maximal likelihood
        ind_good = (lh_array != -np.inf)
        if len(lh_array[ind_good]) == 0:
            return 0, 0
        best = sorted(np.unravel_index(np.argmax(lh_array), lh_array.shape))
        prediction = (int(best[0]), int(best[1]))

        return prediction


def generate_models(min_rep: int, max_rep: int, multiple_bckgs: bool = True) -> Iterator[tuple[int | str, int | str]]:
    """
    Generate all pairs of alleles (models for generation of reads).
    :param min_rep: int - minimal number of repetitions
    :param max_rep: int - maximal number of repetitions
    :param multiple_backgrounds: bool - whether to generate all background states
    :return: generator of allele pairs (numbers or 'E' or 'B')
    """
    for model_index1 in range(min_rep, max_rep):
        for model_index2 in range(model_index1, max_rep):
            yield model_index1, model_index2
        yield model_index1, 'E'
        if multiple_bckgs:
            yield 'B', model_index1

    yield 'B', 'B'
    yield 'E', 'E'


def generate_models_one_allele(min_rep: int, max_rep: int) -> Iterator[tuple[int | str, int | str]]:
    """
    Generate all pairs of alleles (models for generation of reads).
    :param min_rep: int - minimal number of repetitions
    :param max_rep: int - maximal number of repetitions
    :return: generator of allele pairs (numbers or 'E' or 'B'), 'X' for non-existing allele
    """
    for model_index1 in range(min_rep, max_rep):
        yield model_index1, 'X'

    yield 'B', 'X'
    yield 'E', 'X'


def convert_to_sym(max_rep, best: tuple[int, int], monoallelic: bool) -> tuple[int | str, int | str]:
    """
    Convert numeric alleles to their symbolic representations.
    :param best: (int, int) - numeric representation of alleles
    :param monoallelic: bool - if this is monoallelic version
    :return: (int|str, int|str) - symbolic representation of alleles
    """
    # convert it to symbols
    if best[0] == 0 and best[1] == max_rep:
        best_sym = ('E', 'E')
    else:
        def fn1(x):
            return 'E' if x == max_rep else 'B' if x == 0 else x
        best_sym = tuple(map(fn1, best))

    # if mono-allelic return 'X' as second allele symbol
    if monoallelic:
        best_sym = (best_sym[0], 'X')

    return best_sym


def get_confidence(
    lh_array: np.ndarray, predicted: tuple[int, int], max_rep: int, monoallelic: bool = False
) -> Confidences:
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

# DEFAULT_MODEL_PARAMS = (0.001, 0.000105087, 0.0210812, 0.001)
# p1, p2, p3, q = DEFAULT_MODEL_PARAMS
# inserts = q
# deletes = p1 + p2 * n
DEFAULT_MODEL_PARAMS = (0.0001, 0.0001, 0.0, 0.0001)


def model_full(rng: int, n: int) -> np.ndarray:
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

    # print(rng, n)
    p1, p2, _, q = DEFAULT_MODEL_PARAMS
    deletes: np.ndarray = binom.pmf(np.arange(rng), n, clip(1 - (p1 + p2 * n), 0.0, 1.0))
    # print(deletes)
    inserts: np.ndarray = binom.pmf(np.arange(rng), n, q)
    # print(inserts)
    result = combine_distribs(deletes, inserts)
    # print(result)
    return result


def model_bckg(min_rep, max_with_e) -> np.ndarray:
    """Returns ndarray with length ???"""
    result = np.concatenate([
        np.zeros(min_rep, dtype=float),
        np.ones(max_with_e - min_rep, dtype=float) / float(max_with_e - min_rep)
    ])
    return result
