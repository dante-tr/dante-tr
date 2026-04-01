from __future__ import annotations

from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter, ArgumentTypeError
from datetime import datetime
from typing import TypeAlias, Any
from collections import Counter
# from pprint import pprint

import os
import re
import json
import textwrap

import numpy as np
import pandas as pd

from src_new.constants import \
    VERSION, DANTE_DESCRIPTION, MAX_REPETITIONS, \
    MOTIF_COLUMN_NAME, MOTIF_COLUMN_ID, MOTIF_COLUMN_N_MODS, MOTIF_COLUMN_MOD_CLASS, \
    MOTIF_COLUMN_MODULES, MOTIF_COLUMN_MISMATCHES_STR, MOTIF_COLUMN_MODULE_REPETITIONS, \
    MOTIF_COLUMN_MODULE_NOMENCLATURES

from src_new.motif_stats import generate_motif_stats
from src_new.haplotypes import generate_haplotypes


def main() -> None:
    start_time = datetime.now()
    args = load_arguments()

    print(f'DANTE_remaSTR Starting : {start_time:%Y-%m-%d %H:%M:%S}')
    sample = os.path.basename(os.path.normpath(args.output_dir))

    data_json = analyse_motif(args.input_tsv, args.male, sample)

    # os.makedirs(args.output_dir, exist_ok=True)
    # json_example = json.dumps(data_json, indent=2)  # don't sort keys!
    # with open(f"{args.output_dir}/data_v3.json", "w") as f:
    #     f.write(json_example)

    data_json_old = convert_to_old_json(data_json)
    json_example = json.dumps(data_json_old, indent=2)  # don't sort keys!
    with open(f"{args.output_dir}/data.json", "w") as f:
        f.write(json_example)

    end_time = datetime.now()
    print(f'DANTE_remaSTR Stopping : {end_time:%Y-%m-%d %H:%M:%S}')
    print(f'Total time of run      : {end_time - start_time}')


def convert_to_old_json(json):
    return json


def analyse_motif(input_tsv: str, male: bool, sample: str):
    pf = PostFilter()
    main_json: dict[str, Any] = {}
    main_json["dante_version"] = VERSION
    main_json["postfilter_params"] = pf.get_params()
    main_json["sample"] = sample
    main_json["motifs"] = generate_motifs(input_tsv, male)
    return main_json


def generate_motifs(input_tsv: str, male: bool) -> list:
    motif_table = pd.read_csv(input_tsv, sep='\t')

    motif_json: dict[str, Any] = {}
    motif_json["motif_id"] = motif_table[MOTIF_COLUMN_ID].iloc[0]
    motif_json["motif_stats"] = generate_motif_stats(motif_table, male)
    motif_json["nomenclatures"] = generate_nomenclatures(motif_table, male)
    motif_json["modules"] = generate_modules(motif_table, input_tsv, male)
    motif_json["phasings"] = generate_haplotypes(motif_table, male, motif_json["modules"])
    # print(motif_json["phasings"])
    motif_json["phased_seqs"] = generate_phased_seqs(motif_table, male, motif_json["modules"], motif_json["phasings"])
    motif_json["phased_seqs_read_counts"] = generate_read_counts(
        motif_json["nomenclatures"], motif_json["modules"], motif_json["phased_seqs"]
    )

    return [motif_json]


def generate_read_counts(full_nomenclatures: dict, modules: list, phased_seqs: dict) -> dict[str, Any]:
    result: dict[str, Any] = {}
    # get number of reads spanning the whole motif
    h1_count = 0
    h2_count = 0
    for nom in full_nomenclatures:
        if nom["noms"] == phased_seqs["nomenclature1"]:
            h1_count = nom["count"]

    for nom in full_nomenclatures:
        if nom["noms"] == phased_seqs["nomenclature2"]:
            h2_count = nom["count"]
    result["full_reads"] = ((h1_count, phased_seqs["nomenclature1"]), (h2_count, phased_seqs["nomenclature2"]))

    # get number of reads spanning per each module
    module_reads = []
    for i, mod in enumerate(modules):
        module_nomenclatures = mod["nomenclatures"]
        h1_count = 0
        h2_count = 0
        h1_struct = [phased_seqs["nomenclature1"][i]]
        h2_struct = [phased_seqs["nomenclature2"][i]]
        for nom in module_nomenclatures:
            if nom["noms"] == h1_struct:
                h1_count = nom["count"]

        for nom in module_nomenclatures:
            if nom["noms"] == h2_struct:
                h2_count = nom["count"]
        module_reads.append(((h1_count, h1_struct), (h2_count, h2_struct)))
    result["module_reads"] = module_reads

    return result


def generate_nomenclatures(motif_table: pd.DataFrame, male: bool) -> list[dict]:
    motif = Motif(motif_table[MOTIF_COLUMN_NAME].iloc[0], motif_table[MOTIF_COLUMN_ID].iloc[0], male)
    n_modules: int = int(motif_table[MOTIF_COLUMN_N_MODS].iloc[0])
    selected = pd.Series([True] * len(motif_table))
    for module_number in range(1, n_modules - 1):
        x = (motif_table[MOTIF_COLUMN_MOD_CLASS].apply(lambda x: x.split(",")[module_number]) == "Spanning")
        selected &= x

    annotations: list[Annotation] = []
    for _, mt_row in motif_table[selected].iterrows():
        annotations.append(Annotation(mt_row))

    raw_nomenclatures = ["\t".join(annot.module_nomenclatures[1:-1]) for annot in annotations]
    motif_nomenclatures = []
    for (count, location, part) in format_nomenclatures(raw_nomenclatures, motif, None):
        motif_nomenclature: dict[str, Any] = {}
        motif_nomenclature["count"] = count
        motif_nomenclature["location"] = location
        motif_nomenclature["noms"] = part
        motif_nomenclatures.append(motif_nomenclature)
    return motif_nomenclatures


def generate_modules(motif_table: pd.DataFrame, input_tsv: str, male: bool):
    modules = []

    motif = Motif(motif_table[MOTIF_COLUMN_NAME].iloc[0], motif_table[MOTIF_COLUMN_ID].iloc[0], male)
    prediction_json = input_tsv[0:-len(".annotations.tsv")] + ".genotypes.json"
    with open(prediction_json) as f:
        predictions = json.load(f)
    seq = motif.modules_str(include_flanks=True)
    rep_mods = [(int(i), str(seq), int(num)) for i, (seq, num) in enumerate(motif.modules) if num > 1]
    for idx, (module_number, _, _) in enumerate(rep_mods):
        selected = motif_table[MOTIF_COLUMN_MOD_CLASS].apply(lambda x: x.split(",")[module_number])
        anns_spanning = [Annotation(row) for _, row in motif_table[selected == "Spanning"].iterrows()]
        anns_flanking = [Annotation(row) for _, row in motif_table[selected == "Flanking"].iterrows()]
        anns_inrepeat = [Annotation(row) for _, row in motif_table[selected == "In-repeat"].iterrows()]

        heatmap_data, prediction, raw_confidence = do_full_prediction2(motif, anns_spanning, module_number, predictions, idx)

        raw_nomenclature = [annot.module_nomenclatures[module_number] for annot in anns_spanning]
        mod_nomenclatures = format_nomenclatures(raw_nomenclature, motif, None)
        module_nomenclatures = []
        for mod_nomenclature in mod_nomenclatures:
            nomenclature: dict[str, Any] = {}
            nomenclature["count"] = mod_nomenclature[0]
            nomenclature["location"] = mod_nomenclature[1]
            nomenclature["noms"] = mod_nomenclature[2]
            module_nomenclatures.append(nomenclature)

        read_counts = None
        if len(anns_spanning) != 0 or len(anns_flanking) != 0:
            read_counts = write_histogram_image(anns_spanning, anns_flanking, anns_inrepeat, module_number)
        else:
            print(f"Zero reads in {motif.name}")

        graph_data: GraphData
        graph_data = (read_counts, heatmap_data, None)

        # until here it is quite ok, but here the mess begins
        a1_reads: int | str = '---'
        a2_reads: int | str = '---'
        indels: float | str = '---'
        a1_indels: float | str = '---'
        a2_indels: float | str = '---'
        mismatches: float | str = '---'
        a1_mismatches: float | str = '---'
        a2_mismatches: float | str = '---'
        # get info about errors and number of reads from quality annotations if provided
        if anns_spanning is not None:
            # get info about number of reads
            __7__a1 = int(prediction[0]) if isinstance(prediction[0], int) else None
            __7__a2 = int(prediction[1]) if isinstance(prediction[1], int) else None
            a1_reads = 0 if __7__a1 is None else len([a for a in anns_spanning if a.module_repetitions[module_number] == __7__a1])
            a2_reads = 0 if __7__a2 is None else len([a for a in anns_spanning if a.module_repetitions[module_number] == __7__a2])

            # get info about errors
            __7__errors = [a.get_module_errors(motif, module_number) for a in anns_spanning]
            __7__errors_a1 = [a.get_module_errors(motif, module_number) for a in anns_spanning if a.module_repetitions[module_number] == __7__a1]
            __7__errors_a2 = [a.get_module_errors(motif, module_number) for a in anns_spanning if a.module_repetitions[module_number] == __7__a2]
            assert len([l for i, _, l in __7__errors if l == 0]) == 0

            # extract error metrics
            indels, mismatches = errors_per_read(__7__errors, relative=True)
            a1_indels, a1_mismatches = errors_per_read(__7__errors_a1, relative=True)
            a2_indels, a2_mismatches = errors_per_read(__7__errors_a2, relative=True)

        # fill templates:
        confidence = float_to_str(raw_confidence[0], percents=True)
        indels = float_to_str(indels, decimals=2)
        mismatches = float_to_str(mismatches, decimals=2)

        a1_confidence = float_to_str(raw_confidence[1], percents=True)
        a1_indels = float_to_str(a1_indels, decimals=2)
        a1_mismatches = float_to_str(a1_mismatches, decimals=2)

        a2_confidence = float_to_str(raw_confidence[2], percents=True)
        a2_indels = float_to_str(a2_indels, decimals=2)
        a2_mismatches = float_to_str(a2_mismatches, decimals=2)

        myid = list(map(int, str(module_number).split('_')))
        sequence, _ = highlight_subpart(seq, myid)
        __17__motif_name_part1 = f'{motif.name.replace("/", "_")}'
        __17__motif_name_part2 = f'{",".join(map(str, myid)) if myid is not None else "mot"}'
        __17__motif_name_long = f'{__17__motif_name_part1}_{__17__motif_name_part2}'
        __17__motif_clean = re.sub(r'[^\w_]', '', __17__motif_name_long)
        mname = __17__motif_clean.rsplit('_', 1)[0]
        module_id = mname if myid == [1] else __17__motif_clean  # trick to solve static html

        module: dict[str, Any] = {}
        module["module_id"] = module_id
        module["id"] = myid
        module["sequence"] = sequence
        module["nomenclatures"] = module_nomenclatures
        module["allele_1"] = (prediction[0], a1_confidence, a1_indels, a1_mismatches, a1_reads)  # this could contain seq_pred, seq_reads
        module["allele_2"] = (prediction[1], a2_confidence, a2_indels, a2_mismatches, a2_reads)  # -||-
        module["stats"] = (confidence, indels, mismatches)
        module["raw_conf"] = raw_confidence
        module["reads_spanning"] = len(anns_spanning)
        module["reads_flanking"] = len(anns_flanking)
        module["graph_data"] = graph_data
        modules.append(module)
    return modules


def do_full_prediction2(
    motif, anns_spanning, module_number, predictions, idx: int
) -> tuple[None | ProbHeatmap, tuple[int | str, int | str], tuple[float, ...]]:
    module = predictions["modules"][idx]
    spanning_observed_counts = [ann.module_repetitions[module_number] for ann in anns_spanning]

    prediction: tuple[int | str, int | str]
    if len(spanning_observed_counts) == 0:
        likelihoods, prediction, raw_confidence = None, ('B', 'B'), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    else:
        likelihoods = parse_likelihoods(module["likelihoods"])
        prediction = parse_prediction_from_enum(module["predictions_enum"])
        raw_confidence = module["confidences"]

    heatmap_data = None
    if likelihoods is not None:
        max_spanning_reps = max(spanning_observed_counts)
        _max_rep = max_spanning_reps + 1
        _min_rep = get_min_rep(spanning_observed_counts)
        exp_idx = likelihoods.shape[0] - 2
        bkg_idx = likelihoods.shape[0] - 1
        likelihoods = transform_to_old_format(likelihoods, _min_rep, _max_rep, exp_idx, bkg_idx)

        my_min_rep = 1
        while np.isinf(likelihoods[my_min_rep][my_min_rep]):
            my_min_rep += 1
        my_max_rep = likelihoods.shape[0]
        my_max_with_e = likelihoods.shape[1]
        heatmap_data = draw_pcolor(likelihoods, motif.nomenclature, my_min_rep, my_max_rep, my_max_with_e)
    else:
        print(f"Likelihood array is None for {motif.name}.")

    return heatmap_data, prediction, raw_confidence


def parse_prediction_from_enum(json_p_enum) -> tuple[int | str, int | str]:
    if isinstance(json_p_enum[0], dict):
        a = json_p_enum[0]["Num"]
    else:  # it is string
        a = 'E' if json_p_enum[0] == "Expansion" else 'B'

    if isinstance(json_p_enum[1], dict):
        b = json_p_enum[1]["Num"]
    else:  # it is string
        b = 'E' if json_p_enum[1] == "Expansion" else 'B'

    return (a, b)


def parse_likelihoods(json_likelihoods) -> np.ndarray:
    x = np.array(json_likelihoods["data"])
    y: np.ndarray = x.reshape(json_likelihoods["dim"]).astype(float)
    y[np.isnan(y)] = -np.inf
    return y


def get_min_rep(spanning_obs_counts: list[int]) -> int:
    MIN_REPETITIONS = 1
    OVERHEAD = 3
    return max(MIN_REPETITIONS, min(spanning_obs_counts) - OVERHEAD)  # inclusive


# TODO: split this into somethings integratable to class and conversion to old
def transform_to_old_format(lhoods, min_rep, max_rep, exp_idx, bkg_idx):
    # print(lhoods.shape)
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
    # print(likelihoods.shape)
    return likelihoods


def generate_phased_seqs(motif_table: pd.DataFrame, male: bool, modules: dict, phasings: dict):
    motif = Motif(motif_table[MOTIF_COLUMN_NAME].iloc[0], motif_table[MOTIF_COLUMN_ID].iloc[0], male)
    assignment_factor = 0.8
    hap1 = []
    hap2 = []
    nomenclatures = []
    for mod in modules:
        mod_id = mod["id"][0]
        prediction = (mod["allele_1"][0], mod["allele_2"][0])

        selected = motif_table[MOTIF_COLUMN_MOD_CLASS].apply(lambda x: x.split(",")[mod_id])
        spanning = [Annotation(row) for _, row in motif_table[selected == "Spanning"].iterrows()]
        hap1.append(str(prediction[0]))
        hap2.append(str(prediction[1]))

        raw_nomenclatures = [annot.module_nomenclatures[mod_id] for annot in spanning]
        nomenclature = list(map(nom_count_to_triple, format_nomenclatures(raw_nomenclatures, motif, None)))
        nomenclatures.append(nomenclature)

    phase1 = []
    phase2 = []
    for phasing in phasings:
        phase1.append(phasing["allele_1"][0])
        phase2.append(phasing["allele_2"][0])

    hap1_full, hap2_full, err1, err2 = phase_full_locus(hap1, hap2, phase1, phase2)
    nomenclature1, nomenclature1_len, nomenclatures, errs1 = augment_nomenclature(motif, hap1_full, nomenclatures, assignment_factor)
    nomenclature2, nomenclature2_len, nomenclatures, errs2 = augment_nomenclature(motif, hap2_full, nomenclatures, assignment_factor)
    result: dict[str, Any] = {}
    result["motif_name"] = motif.name
    result["nomenclature1"] = nomenclature1
    result["nomenclature1_occ"] = hap1_full
    result["nomenclature1_len"] = nomenclature1_len
    result["errors1"] = err1 + errs1
    result["nomenclature2"] = nomenclature2
    result["nomenclature2_occ"] = hap2_full
    result["nomenclature2_len"] = nomenclature2_len
    result["errors2"] = err2 + errs2
    return result


def augment_nomenclature(
    motif: Motif, hapl: list[str], nomenclatures: list[list[tuple[int, int, str]]], assignment_factor: float
) -> tuple[list[str], int, list[list[tuple[int, int, str]]], list[str]]:

    errors = set()
    result = []
    for i, x in enumerate(hapl):
        try:
            count = int(x)
        except ValueError:
            errors.add("contains B/E/X")
            result.append(f"err[{x}]")
            continue

        assigned = False
        for j, nom in enumerate(nomenclatures[i]):
            num_repeats, occ, representation = nom
            mod_seq = motif.modules[1 + i][0]
            rep_len = hgvs_to_len(representation)
            mod_len = len(mod_seq) * count

            if num_repeats == count or rep_len == mod_len:
                result.append(representation)
                occ = int(occ * (1.0 - assignment_factor))
                nomenclatures[i][j] = (num_repeats, occ, representation)
                assigned = True
                nomenclatures[i] = sorted(nomenclatures[i], key=lambda x: -x[1])
                break

        if not assigned:
            errors.add("nomenclature mismatch")
            result.append(f"err[{count}]")

    aug_nom = motif.augmented_nomenclature(result)
    nom_len = sum((hgvs_to_len(aug_nom[i]) for i in range(len(aug_nom))))
    return aug_nom, nom_len, nomenclatures, list(errors)


def hgvs_to_len(hgvs: str) -> int:
    "Used for sorting based on nomenclature length."
    length = 0
    for seq, num in re.findall(r'([A-Z]+)\[([0-9BE]+)\]', hgvs):
        if num == "B":
            length += -1000
        elif num == "E":
            length += 1000
        else:
            length += len(seq) * int(num)

    return length


def nom_count_to_triple(nomenclature: tuple[int, str, list[str]]) -> tuple[int, int, str]:
    count, _, repr_list = nomenclature
    assert len(repr_list) == 1, "This function converts only one representation"
    representation = repr_list[0]
    length = sum(int(num) for _, num in re.findall(r'([A-Z]+)\[(\d+)', representation))
    return length, count, representation


def phase_full_locus(
    h1_full: list[str], h2_full: list[str], hp1: list[str], hp2: list[str]
) -> tuple[list[str], list[str], list[str], list[str]]:
    # ['9', '15'] ['12', '16'] ['9|16'] ['12|15'] -> ['9', '16'] ['12', '15']

    n_diff: int = sum(map(lambda x: x[0] != x[1], zip(h1_full, h2_full)))
    if n_diff <= 1:
        return (h1_full, h2_full, [], [])  # there is nothing to phase

    errors1 = set()
    errors2 = set()
    # print()
    # print(h1_full, h2_full, hp1, hp2)
    different_prefix = False
    for i in range(len(h1_full) - 1):
        if h1_full[i] == h2_full[i]:
            if different_prefix:
                errors1.add("homozygous link")
                errors2.add("homozygous link")
                print(f"Warning: Cannot phase prefix and suffix at position {i}.\n{h1_full} {h2_full}\n{hp1} {hp2}")
            continue
        different_prefix = True

        h1_cur, _ = hp1[i].split("|")
        if h1_cur != h1_full[i]:
            hp1[i], hp2[i] = hp2[i], hp1[i]

        h1_cur, h1_next = hp1[i].split("|")
        _, h2_next = hp2[i].split("|")

        if h1_full[i] == h1_cur and h1_full[i + 1] == h1_next:
            # if h2_full[i] != h2_cur or h2_full[i + 1] != h2_next: WARN?
            pass
        elif h1_full[i] == h1_cur and h1_full[i + 1] == h2_next:
            # if h2_full[i] != h2_cur or h2_full[i + 1] != h1_next: WARN?
            h1_full[i + 1], h2_full[i + 1] = h2_full[i + 1], h1_full[i + 1]
        else:
            print(f"Warning: {h1_full[i:i + 2]} {h2_full[i:i + 2]} {hp1[i]} {hp2[i]} is inconsistent.")

    for i in range(len(h1_full) - 1):
        new_p1 = f"{h1_full[i]}|{h1_full[i + 1]}"
        new_p2 = f"{h2_full[i]}|{h2_full[i + 1]}"
        from_genotyping = sorted([new_p1, new_p2])
        from_phasing = sorted([hp1[i], hp2[i]])
        if from_genotyping != from_phasing:
            errors1.add("genotyping-phasing inconsistence")
            errors2.add("genotyping-phasing inconsistence")

    # print(f"-> {h1_full} {h2_full}\n")
    return (h1_full, h2_full, list(errors1), list(errors2))


def format_nomenclatures(
    nomenclatures: list[str], motif: Motif, nomenclature_limit: int | None = None
) -> list[tuple[int, str, list[str]]]:
    counter = sorted(Counter(nomenclatures).items(), key=lambda k: (-k[1], k[0]))
    lines = []
    for nomenclature, count in counter:
        ref = f'{motif.chrom}:g.{motif.start}_{motif.end}'
        parts = nomenclature.rstrip().split('\t')

        lines.append((count, ref, parts))

        if nomenclature_limit is not None and len(lines) >= nomenclature_limit:
            break

    return lines


def load_arguments() -> Namespace:
    """
    Loads and parses the arguments.
    :return: args - parsed arguments
    """
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=textwrap.dedent(DANTE_DESCRIPTION)
    )

    options = parser.add_argument_group('Options')
    options.add_argument(
        '--input-tsv', '-i', default="",
        help='Input annotation table as obtained by remaSTR. Default=stdin'
    )
    options.add_argument(
        '--output-dir', '-o', type=str, default="dante_out",
        help='Output destination (directory). Default=./dante_out/'
    )
    options.add_argument(
        '--verbose', '-v', action='store_true',
        help='Print all the outputs. Default is to print only the result table to stdout.'
    )
    options.add_argument(
        '--male', action='store_true',
        help='Indicate that the sample is male. Process motifs from chrX/chrY as mono-allelic.'
    )
    options.add_argument(
        '--nomenclatures', '-n', type=positive_int, default=5,
        help='Number of nomenclature strings to add to reports. Default=5'
    )
    options.add_argument(
        '--cutoff-alignments', type=positive_int, default=20,
        help='How many bases to keep beyond annotated part. Default=20'
    )

    args = parser.parse_args()

    return args


def positive_int(value: str) -> int:
    try:
        int_value = int(value)
    except ValueError:
        raise ArgumentTypeError(f'Value {value} is not integer') from None
    if int_value < 0:
        raise ArgumentTypeError(f'Value {value} is negative')
    return int_value


class Motif:
    """
    Class to represent DNA motifs.

    :ivar chrom: Chromosome name.
    :ivar start: Start position of the motif.
    :ivar end: End position of the motif.
    :ivar modules: A list of tuples containing sequence and repetition count.
    :ivar name: name of the Motif
    :ivar motif: motif nomenclature
    """

    def __init__(self, motif: str, name: str | None, male: bool) -> None:
        """
        Initialize a Motif object.
        :param motif: The motif string in the format "chrom:start_end[A][B]..."
        :param name: optional name of the motif
        """
        # remove whitespace
        nomenclature = motif.strip().replace(' ', '')
        name = (name if name is not None else nomenclature).replace(':', '-').replace('.', '_').replace('/', '_')

        # extract prefix, first number, second number
        tmp = re.match(r'([^:]+):g\.(\d+)_(\d+)(.*)', nomenclature)
        if tmp is None:
            raise ValueError(f"{nomenclature} has incorrect format")
        chrom, start, end, remainder = tmp.groups()

        # extract sequence and repetition count
        modules = [(str(seq), int(num)) for seq, num in re.findall(r'([A-Z]+)\[(\d+)', remainder)]
        modules = [('left_flank', 1)] + modules + [('right_flank', 1)]

        # store members
        self.nomenclature: str = nomenclature
        self.name: str = name
        self.chrom: str = chrom
        self.start: int = int(start)
        self.end: int = int(end)
        self.modules: list[tuple[str, int]] = modules

    def __getitem__(self, index: int) -> tuple[str, int]:
        """
        Returns module at given index.
        :param index: The index of the module to fetch.
        :return: The module at the given index.
        """
        return self.modules[index]

    def __str__(self) -> str:
        """
        Returns string representation of the Motif object.
        :return: String representation in the format "chrom:start_end[A][B]..."
        """
        return f'{self.chrom}:g.{self.start}_{self.end}' + self.modules_str(include_flanks=False)

    def __lt__(self, obj):
        """
        Less than for sorting purposes.
        :return: bool - if this object comes before the other
        """
        return self.name < obj.name

    def __eq__(self, obj):
        """
        Equal to for sorting purposes.
        :return: bool - if this object is equal to the other
        """
        return self.name == obj.name

    def augmented_nomenclature(self, rep_counts: list[str]) -> list[str]:
        modules = []
        i = 0
        for seq, num in self.modules[1:-1]:
            if num == 1:
                modules.append(f"{seq}[{num}]")
            else:
                if rep_counts[i].startswith("err"):
                    x = rep_counts[i][4:-1]
                    modules.append(f"{seq}[{x}]")
                else:
                    modules.append(rep_counts[i])
                i += 1
        assert i == len(rep_counts), "Invalid augmentation."
        return modules

    def modules_str(self, include_flanks: bool = False) -> str:
        """
        Returns string representation of modules
        :param include_flanks: bool - include flank modules?
        :return: String representation of modules
        """
        if include_flanks:
            return ''.join([f'{seq}[{num}]' for seq, num in self.modules])
        return ''.join([f'{seq}[{num}]' for seq, num in self.modules[1:-1]])


class Annotation:
    def __init__(self, row: pd.Series):
        # Store arguments into instance variables
        self.states = row[MOTIF_COLUMN_MODULES]

        # Calculate insertion/deletion/mismatch string
        self.mismatches_string = row[MOTIF_COLUMN_MISMATCHES_STR]

        # Calculate number of insertions, deletions and normal bases

        # Number of STR motif repetitions and sequences of modules
        self.module_repetitions = list(map(int, row[MOTIF_COLUMN_MODULE_REPETITIONS].split(",")))
        self.module_nomenclatures = list(row[MOTIF_COLUMN_MODULE_NOMENCLATURES].split(","))

    def get_module_errors(self, motif: Motif, module_num: int, overhang: int | None = None) -> tuple[int, int, int]:
        # get overhang as module length
        if overhang is None:
            seq, _ = motif.modules[module_num]
            overhang = len(seq)

        # define module character
        char_to_search = chr(ord('0') + module_num)

        # if the annotation does not have this module, return 0
        if char_to_search not in self.states:
            return 0, 0, 0

        # search for the annotation of the module
        start = max(0, self.states.find(char_to_search) - overhang)
        end = min(self.states.rfind(char_to_search) + overhang + 1, len(self.states))

        # count errors
        indels = self.mismatches_string[start:end].count('I') + self.mismatches_string[start:end].count('D')
        mismatches = self.mismatches_string[start:end].count('M')

        # return indels, mismatches, and length
        return indels, mismatches, end - start

    def primers(self, index_rep: int) -> int:
        """
        Count how any primers it has on repetition index.
        :param index_rep: int - index of the repetition, that we are looking at
        :return: int - number of primers (0-2)
        """
        primers = 0
        if index_rep > 0 and self.module_repetitions[index_rep - 1] > 0:
            primers += 1
        if index_rep + 1 < len(self.module_repetitions) and self.module_repetitions[index_rep + 1] > 0:
            primers += 1
        return primers


def errors_per_read(
    errors: list[tuple[int, int, int]], relative: bool = False
) -> tuple[float | str, float | str]:
    """
    Count number of errors per read. Relative per length or absolute number.
    :param errors: list[tuple[int, int, int]] - indels, mismatches and length of module
    :param relative: bool - relative?
    :return: tuple[float, float] - number of indels, mismatches per hundred reads
    """
    # if we have no reads, return '---'
    if len(errors) == 0:
        return '---', '---'

    if relative:
        mean_length = np.mean([length for _, _, length in errors])
        return (
            float(np.mean([indels / float(length) for indels, _, length in errors]) * mean_length),
            float(np.mean([mismatches / float(length) for _, mismatches, length in errors]) * mean_length)
        )
    return (
        float(np.mean([indels for indels, _, _ in errors])),
        float(np.mean([mismatches for _, mismatches, _ in errors]))
    )


class PostFilter:
    MIN_REP_LEN = 3
    MIN_REP_CNT = 1
    MAX_ABS_ERROR = None
    MAX_REL_ERROR = 1.0

    def __init__(self):
        self.min_rep_len = self.MIN_REP_LEN
        self.min_rep_cnt = self.MIN_REP_CNT
        self.max_rel_error = self.MAX_REL_ERROR
        self.max_abs_error = self.MAX_ABS_ERROR

    def get_params(self):
        pf_error = f'{self.max_rel_error * 100:.0f}%'
        if self.max_abs_error is not None:
            pf_error += f' (abs={self.max_abs_error})'
        return self.min_rep_len, self.min_rep_cnt, pf_error


def sorted_repetitions(annotations: list[Annotation]) -> list[tuple[tuple[int, ...], int]]:
    """
    Aggregate same repetition counts for annotations and sort them according to quantity of repetitions of each module
    :param annotations: Annotated reads
    :return: list of (repetitions, count), sorted by repetitions
    """
    count_dict = Counter(tuple(annot.module_repetitions) for annot in annotations)
    return sorted(count_dict.items(), key=lambda k: k[0])


def write_histogram_image(
        annotations: list[Annotation], filt_annot: list[Annotation], anns_inrepeat: list[Annotation], index_rep: int
) -> HistReadCounts:
    """
    Stores quantity of different combinations of module repetitions, generates separate graph image for each module
    :param out_prefix: Output file prefix
    :param annotations: Annotated reads.
    :param filt_annot: Annotated reads (filtered)
    :param index_rep: int - index of repetition module of a motif
    """
    repetitions = sorted_repetitions(annotations)
    repetitions_filt = sorted_repetitions(filt_annot)
    repetitions_inrep = sorted_repetitions(anns_inrepeat)

    spanning_counts = [(r[index_rep], c) for r, c in repetitions]
    filtered_counts = [(r[index_rep], c) for r, c in repetitions_filt]
    inread_counts = [(r[index_rep], c) for r, c in repetitions_inrep]

    xm = max(
        [r for r, c in spanning_counts]
        + [r for r, c in filtered_counts]
        + [r for r, c in inread_counts]
        + [MAX_REPETITIONS]
    )

    # set data
    spanning = [0] * (xm + 1)
    for r, c in spanning_counts:
        spanning[r] += c

    flanking = spanning.copy()
    for r, c in filtered_counts:
        flanking[r] += c

    inread = flanking.copy()
    for r, c in inread_counts:
        inread[r] += c

    # plot_histogram_image_plotly(out_prefix, spanning, flanking, inread)

    only_flanking = [df - d for df, d in zip(flanking, spanning)]
    only_inread = [di - df for di, df in zip(inread, flanking)]
    return spanning, only_flanking, only_inread


def draw_pcolor(
    lh_array: np.ndarray, name: str, min_rep: int, max_rep: int, max_with_e: int, lognorm: bool = True
) -> ProbHeatmap:

    ind_good = (lh_array < 0.0) & (lh_array > -1e10) & (lh_array != np.nan)
    z_min, z_max = min(lh_array[ind_good]), max(lh_array[ind_good])
    max_str = len(lh_array)
    if lognorm:
        lh_view = -np.log(-lh_array)
        z_min = -np.log(-z_min)
        z_max = -np.log(-z_max)
    else:
        lh_view = lh_array.copy()

    # background (B, i) - copy it below min_rep
    lh_view[min_rep - 1, :] = lh_view[0, :]

    lh_copy = lh_view.copy()
    lh_copy[-1, min_rep] = lh_copy[0, 0]
    lh_copy[-1, min_rep + 1] = lh_copy[0, max_rep]

    title = '%s likelihood of options (%s)' % ('Loglog' if lognorm else 'Log', name)
    # print(lh_copy.shape, lognorm, title, max_str)
    heatmap_data = save_pcolor_plotly_file(lh_copy, lognorm, title, max_str, min_rep, max_with_e)
    return heatmap_data


def save_pcolor_plotly_file(
    lh_copy: np.ndarray, lognorm: bool, title: str, max_str: int, min_rep: int, max_with_e: int,
    start_ticks: int = 5, step_ticks: int = 5
) -> ProbHeatmap:
    text = [['' for _ in range(max_str - min_rep + 1)] for _ in range(max_str - min_rep + 1)]
    text[-1][0] = 'B'
    text[-1][1] = 'E'

    hovertext = []
    for j in ['B'] + list(range(min_rep, max_str)):
        inner = [f'{j}/{i}' for i in list(range(min_rep, max_str)) + ['E']]
        hovertext.append(inner)

    hovertext[0][-1] = 'E/E'
    hovertext[-1][0] = 'B'
    hovertext[-1][1] = 'E'

    z: list[list[float | None]] = lh_copy[min_rep - 1:, min_rep:].tolist()
    for i in range(len(z)):
        for j in range(len(z[0])):
            if z[i][j] == -np.inf:
                z[i][j] = None
    hovertext2: list[list[str]] = hovertext
    y_tickvals: list[int] = list(range(start_ticks - min_rep + 1, max_str - min_rep + 1, step_ticks)) + [0]
    y_ticktext: list[int | str] = list(range(start_ticks, max_str, step_ticks)) + ['B']
    x_tickvals: list[int] = list(range(start_ticks - min_rep, max_str - min_rep, step_ticks)) + [int(max_str - min_rep)]
    x_ticktext: list[int | str] = list(range(start_ticks, max_str, step_ticks)) + ['E(>%d)' % (max_with_e - 2)]
    x_pos: float = max_str - min_rep - 0.5

    return (z, hovertext2, y_tickvals, y_ticktext, x_tickvals, x_ticktext, x_pos)


def highlight_subpart(seq: str, highlight: int | list[int]) -> tuple[str, str]:
    """
    Highlights subpart of a motif sequence
    :param seq: str - motif sequence
    :param highlight: int/list(int) - part ot highlight
    :return: str, str - motif sequence with highlighted subpart, highlighted subpart
    """
    if highlight is None:
        return seq, ''

    str_part = []
    highlight1 = np.array(highlight)
    split = [f'{s}]' for s in seq.split(']') if s != '']
    for h in highlight1:
        str_part.append(split[h])
        split[h] = f'<b><u>{split[h]}</u></b>'
    return ''.join(split), ''.join(str_part)


def float_to_str(c: float | str, percents: bool = False, decimals: int = 1) -> str:
    """
    Convert float confidence to string.
    :param c: float/str - confidence
    :param percents: bool - whether to output as a percents or not
    :param decimals: int - how many decimals to round to
    :return: str - converted to string
    """
    if isinstance(c, float):
        return f'{c * 100: .{decimals}f}%' if percents else f'{c: .{decimals}f}'
    return c


Hist2DGraph: TypeAlias = tuple[list[list[int]], list[list[int]], list[list[str]], str, str]
HistReadCounts: TypeAlias = tuple[list[int], list[int], list[int]]
ProbHeatmap: TypeAlias = tuple[list[list[float | None]], list[list[str]], list[int], list[int | str], list[int], list[int | str], float]
GraphData: TypeAlias = tuple[HistReadCounts | None, ProbHeatmap | None, Hist2DGraph | None]


# %%
if __name__ == '__main__':
    main()
