from src_new.genotyping import genotype as genotype1
from src_new.genotyping2 import genotype as genotype2
import numpy as np


def test_biallelic_genotype_function() -> None:
    spanning_obsed_cs = [15, 15, 15, 15, 16, 15, 16, 15, 16, 16, 16, 15, 15, 15, 15, 15, 15, 15, 16, 16, 15, 15, 16, 15, 15, 16, 16, 15, 16, 15, 15, 16, 15, 15, 16, 15, 15, 16, 16, 15, 16, 15, 15, 8, 7, 7, 7, 7]
    spanning_read_len = [148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148]
    flanking_obsed_cs = [1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 10, 10, 10, 11, 11, 12, 12, 12, 14, 14, 14, 14, 14, 14, 16, 16, 16, 17, 14, 16, 16, 14, 14, 14, 14, 12, 11, 11, 10, 8, 10, 9, 9, 8, 7, 7, 6, 6, 5, 4, 4, 3, 3, 2, 2, 2, 2]
    flanking_read_len = [148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148]
    monoallelic_motif = False

    likelihood_matrix1, prediction1, confidence1 = genotype1(spanning_obsed_cs, spanning_read_len, flanking_obsed_cs, flanking_read_len, monoallelic_motif, 1, 3, 3)
    likelihood_matrix2, prediction2, confidence2 = genotype2(spanning_obsed_cs, spanning_read_len, flanking_obsed_cs, flanking_read_len, monoallelic_motif, 1, 3, 3)

    # biallelic
    assert prediction1 == prediction2
    # assert np.allclose(confidence1, confidence2, equal_nan=True)
    assert confidence1 == confidence2
    # assert np.allclose(likelihood_matrix1, likelihood_matrix2, equal_nan=True)
    assert (likelihood_matrix1 == likelihood_matrix2).all()


def test_monoallelic_genotype_function() -> None:
    spanning_obsed_cs = [15, 15, 15, 15, 16, 15, 16, 15, 16, 16, 16, 15, 15, 15, 15, 15, 15, 15, 16, 16, 15, 15, 16, 15, 15, 16, 16, 15, 16, 15, 15, 16, 15, 15, 16, 15, 15, 16, 16, 15, 16, 15, 15, 8, 7, 7, 7, 7]
    spanning_read_len = [148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148]
    flanking_obsed_cs = [1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 10, 10, 10, 11, 11, 12, 12, 12, 14, 14, 14, 14, 14, 14, 16, 16, 16, 17, 14, 16, 16, 14, 14, 14, 14, 12, 11, 11, 10, 8, 10, 9, 9, 8, 7, 7, 6, 6, 5, 4, 4, 3, 3, 2, 2, 2, 2]
    flanking_read_len = [148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148]
    monoallelic_motif = True

    likelihood_matrix1, prediction1, confidence1 = genotype1(spanning_obsed_cs, spanning_read_len, flanking_obsed_cs, flanking_read_len, monoallelic_motif, 1, 3, 3)
    likelihood_matrix2, prediction2, confidence2 = genotype2(spanning_obsed_cs, spanning_read_len, flanking_obsed_cs, flanking_read_len, monoallelic_motif, 1, 3, 3)

    # monoallelic
    assert prediction1 == prediction2
    # assert np.allclose(confidence1, confidence2, equal_nan=True)
    assert confidence1 == confidence2
    # assert np.allclose(likelihood_matrix1, likelihood_matrix2, equal_nan=True)
    assert (likelihood_matrix1 == likelihood_matrix2).all()


def test_BSS_motif() -> None:
    spanning_obsed_cs = [9, 12, 12, 12, 12, 12, 12, 12, 1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 4, 2, 3, 1, 2, 1, 1, 1, 4]
    spanning_read_len = [148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 149, 148, 149, 148, 148, 148, 148, 148, 148]
    flanking_obsed_cs = [1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 4, 6, 6, 6, 6, 7, 7, 6, 7, 8, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10, 11, 11, 10, 10, 11, 11, 11, 11, 11, 10, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 16, 15, 16, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21, 22, 23, 24, 24, 25, 24, 2, 6, 2, 1, 1, 13, 2, 39, 1, 1, 1, 36, 2, 2, 16]
    flanking_read_len = [148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 149, 148, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 148, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 148, 148, 148, 149, 148, 149, 148, 149, 148, 148, 148, 149, 148, 148, 149]
    monoallelic_motif = False
    # 1 16 41

    likelihood_matrix1, prediction1, confidence1 = genotype1(spanning_obsed_cs, spanning_read_len, flanking_obsed_cs, flanking_read_len, monoallelic_motif, 1, 3, 3)
    likelihood_matrix2, prediction2, confidence2 = genotype2(spanning_obsed_cs, spanning_read_len, flanking_obsed_cs, flanking_read_len, monoallelic_motif, 1, 3, 3)

    # BSS
    print(prediction1, prediction2)
    assert prediction1 == prediction2
    # assert np.allclose(confidence1, confidence2, equal_nan=True)
    assert confidence1 == confidence2
    # assert np.allclose(likelihood_matrix1, likelihood_matrix2, equal_nan=True)
    assert (likelihood_matrix1 == likelihood_matrix2).all()

# end
