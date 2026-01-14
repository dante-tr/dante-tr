from src_new.genotyping import genotype
from src_new.genotyping2 import genotype as genotype2
import numpy as np


def test_genotype_function() -> None:
    spanning_observed_counts = [
        15, 15, 15, 15, 16, 15, 16, 15, 16, 16, 16, 15, 15, 15, 15, 15, 15, 15,
        16, 16, 15, 15, 16, 15, 15, 16, 16, 15, 16, 15, 15, 16, 15, 15, 16, 15,
        15, 16, 16, 15, 16, 15, 15, 8, 7, 7, 7, 7
    ]
    spanning_read_lengths = [
        148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148,
        148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148,
        148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148,
        148, 148, 148, 148, 148, 148
    ]
    flanking_observed_counts = [
        1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 6, 6,
        6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 10, 10, 10, 11, 11, 12, 12, 12, 14, 14,
        14, 14, 14, 14, 16, 16, 16, 17, 14, 16, 16, 14, 14, 14, 14, 12, 11, 11,
        10, 8, 10, 9, 9, 8, 7, 7, 6, 6, 5, 4, 4, 3, 3, 2, 2, 2, 2
    ]
    flanking_read_lengths = [
        148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148,
        148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148,
        148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148,
        148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148,
        148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148,
        148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148
    ]
    monoallelic_motif = False

    likelihood_matrix, prediction, confidence = genotype(
        spanning_observed_counts, spanning_read_lengths,
        flanking_observed_counts, flanking_read_lengths,
        monoallelic_motif, 1, 3, 3
    )

    likelihood_matrix2, prediction2, confidence2 = genotype2(
        spanning_observed_counts, spanning_read_lengths,
        flanking_observed_counts, flanking_read_lengths,
        monoallelic_motif, 1, 3, 3
    )

    assert prediction == prediction2
    assert confidence == confidence2
    assert (likelihood_matrix == likelihood_matrix2).all()

    # 2026-01-13
    # print(prediction, confidence)
    # np.save('test_genotype_function.npy', likelihood_matrix)
    # (15, 16)
    # (np.float64(1.0), np.float64(1.0), np.float64(1.0),
    #  np.float64(5.7556624067001616e-71), np.float64(1.5837631846830353e-25),
    #  np.float64(7.282452029483091e-179), np.float64(9.043041101350901e-40))
    # cmp test_genotype_function.npy test_genotype_function_bup.npys

