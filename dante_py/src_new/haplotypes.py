from __future__ import annotations

import re
from collections import Counter
from typing import TypeAlias

import pandas as pd
import numpy as np
from src_new.constants import MOTIF_COLUMN_NAME, MOTIF_COLUMN_ID, MOTIF_COLUMN_MOD_CLASS, MAX_REPETITIONS

Hist2DGraph: TypeAlias = tuple[list[list[int]], list[list[int]], list[list[str]], str, str]


def new_phase(
    motif_table: pd.DataFrame, prev_module_num: int, curr_module_num: int
) -> tuple[tuple[str, str], tuple[str, str, str]]:
    mod1 = motif_table[MOTIF_COLUMN_MOD_CLASS].apply(lambda x: x.split(",")[prev_module_num])
    mod2 = motif_table[MOTIF_COLUMN_MOD_CLASS].apply(lambda x: x.split(",")[curr_module_num])
    df_2good = motif_table[(mod1 == "Spanning") & (mod2 == "Spanning")]

    n = len(df_2good)
    if n == 0:
        return ('-|-', '-|-'), ('-/0', '-/0', '-/0')

    reps = []
    for _, r in df_2good.iterrows():
        x = list(map(int, r["module_repetitions"].split(",")))
        reps.append((x[prev_module_num], x[curr_module_num]))

    repetitions = Counter(reps)

    # pick the highest two
    most_common = repetitions.most_common(2)
    rep1, cnt1 = most_common[0]
    rep2, cnt2 = most_common[1] if len(most_common) >= 2 else (('-', '-'), 0)

    # output phasing with number of supported reads
    phasing = (f'{rep1[0]}|{rep1[1]}', f'{rep2[0]}|{rep2[1]}')
    supported_reads = (f'{cnt1 + cnt2}/{n}', f'{cnt1}/{n}', f'{cnt2}/{n}')

    return phasing, supported_reads


def generate_haplotypes(motif_table: pd.DataFrame, male: bool, modules2: dict) -> list:
    motif = Motif(motif_table[MOTIF_COLUMN_NAME].iloc[0], motif_table[MOTIF_COLUMN_ID].iloc[0], male)
    seq = motif.modules_str(include_flanks=True)

    modules = []

    __5__repeating_modules = [(int(i), str(seq), int(num)) for i, (seq, num) in enumerate(motif.modules) if num > 1]
    for __5__i, (curr_module_num, _, _) in enumerate(__5__repeating_modules[1:], start=1):
        prev_module_num = __5__repeating_modules[__5__i - 1][0]
        __5__phasing, __5__supp_reads = new_phase(motif_table, prev_module_num, curr_module_num)

        mod1 = motif_table[MOTIF_COLUMN_MOD_CLASS].apply(lambda x: x.split(",")[prev_module_num])
        mod2 = motif_table[MOTIF_COLUMN_MOD_CLASS].apply(lambda x: x.split(",")[curr_module_num])

        df_2good = motif_table[(mod1 == "Spanning") & (mod2 == "Spanning")]
        good_left = (mod1 == "Spanning") & (mod2 != "Spanning")
        good_rght = (mod1 != "Spanning") & (mod2 == "Spanning")
        df_1good = motif_table[(good_left) | (good_rght)]
        df_0good = motif_table[(mod1 != "Spanning") & (mod2 != "Spanning")]

        anns_2good = [Annotation(r) for _, r in df_2good.iterrows()]
        anns_1good = [Annotation(r) for _, r in df_1good.iterrows()]
        anns_0good = [Annotation(r) for _, r in df_0good.iterrows()]
        __9__phase = (
            curr_module_num, anns_2good, anns_1good, anns_0good, __5__phasing, __5__supp_reads, prev_module_num
        )

        __9__locus_data2 = generate_locus_data2(__9__phase, motif, seq, motif.name, None)

        __1__module = {}
        __1__module["phasing_id"] = __9__locus_data2[0]
        __1__module["ids"] = __9__locus_data2[1]
        __1__module["sequence"] = __9__locus_data2[2]

        module_nomenclatures = []
        for __1__old_nomenclature in __9__locus_data2[3]:
            __1__nomenclature = {}
            __1__nomenclature["count"] = __1__old_nomenclature[0]
            __1__nomenclature["location"] = __1__old_nomenclature[1]
            __1__nomenclature["noms"] = __1__old_nomenclature[2]
            module_nomenclatures.append(__1__nomenclature)
        __1__module["nomenclatures"] = module_nomenclatures

        __1__module["allele_1"] = __9__locus_data2[4]
        __1__module["allele_2"] = __9__locus_data2[5]
        __1__module["stats"] = __9__locus_data2[6]
        __1__module["raw_conf"] = __9__locus_data2[7]
        __1__module["reads_spanning"] = __9__locus_data2[8]
        __1__module["reads_flanking"] = __9__locus_data2[9]
        __1__module["graph_data"] = __9__locus_data2[10]
        modules.append(__1__module)
    return modules


def generate_locus_data2(ph, motif, seq, motif_id, nomenclature_limit):
    if ph is None:
        raise ValueError

    (module_number, anns_2good, anns_1good, anns_0good, phasing1, supp_reads, prev_module_num) = ph
    second_module_number = module_number
    suffix = f'{prev_module_num}_{second_module_number}'

    row = generate_result_line(
        motif, phasing1, supp_reads, len(anns_2good), len(anns_1good), len(anns_0good),
        prev_module_num, second_module_number=module_number
    )
    row_tuple = generate_row(seq, row)

    tmp = generate_motifb64(seq, row)
    (locus_id, highlight, _, _, _, _, sequence) = tmp

    annotations = anns_2good + anns_1good
    if len(annotations) == 0:
        print(f"{motif_id} {suffix} is empty")

    __7__nomenclatures = [
        "\t".join([annot.module_nomenclatures[prev_module_num], annot.module_nomenclatures[second_module_number]])
        for annot in anns_2good
    ]

    locus_nomenclatures = format_nomenclatures(__7__nomenclatures, motif, nomenclature_limit)
    hist2d_data = write_histogram_image2d(
        annotations, prev_module_num, second_module_number,
        motif.module_str(prev_module_num), motif.module_str(second_module_number)
    )

    graph_data = (None, None, hist2d_data)

    _, _, a1_prediction, a1_confidence, a1_reads, a1_indels, a1_mismatches, \
        a2_prediction, a2_confidence, a2_reads, a2_indels, a2_mismatches, \
        confidence, indels, mismatches, spanning_reads, flanking_reads = row_tuple
    raw_confidence = "tmp"
    locus_data2 = (
        locus_id, highlight, sequence, locus_nomenclatures,
        (a1_prediction, a1_confidence, a1_indels, a1_mismatches, a1_reads),
        (a2_prediction, a2_confidence, a2_indels, a2_mismatches, a2_reads),
        (confidence, indels, mismatches),
        raw_confidence,
        spanning_reads, flanking_reads,
        graph_data,
    )
    return locus_data2


def write_histogram_image2d(
    deduplicated: list[Annotation], index_rep: int, index_rep2: int, seq: str, seq2: str
) -> Hist2DGraph | None:
    if deduplicated is None or len(deduplicated) == 0:
        return None

    dedup_reps: list[tuple[tuple[bool, int], tuple[bool, int]]] = []
    for x in deduplicated:
        r_1 = x.get_str_repetitions(index_rep)
        r_2 = x.get_str_repetitions(index_rep2)
        if r_1 is not None and r_2 is not None:
            dedup_reps.append((r_1, r_2))

    if len(dedup_reps) == 0:
        return None

    # assign maximals
    xm = max(r for (_, r), _ in dedup_reps)
    ym = max(r for _, (_, r) in dedup_reps)
    max_ticks = max(ym, xm) + 2
    xm = max(MAX_REPETITIONS, xm)
    ym = max(MAX_REPETITIONS, ym)

    # create data containers
    data = np.zeros((xm + 1, ym + 1), dtype=int)
    data_primer = np.zeros((xm + 1, ym + 1), dtype=int)
    for ((c1, r1), (c2, r2)) in dedup_reps:
        if c1 and c2:
            data[r1, r2] += 1
        if c1 and not c2:
            data_primer[r1, r2:] += 1
        if not c1 and c2:
            data_primer[r1:, r2] += 1

    str1 = 'STR %d [%s]' % (index_rep + 1, seq.split('-')[-1])
    str2 = 'STR %d [%s]' % (index_rep2 + 1, seq2.split('-')[-1])

    def parse_labels(num, num_primer):
        if num == 0 and num_primer == 0:
            return ''
        if num == 0 and num_primer != 0:
            return '0/%s' % str(num_primer)
        if num != 0 and num_primer == 0:
            return '%s/0' % str(num)
        return '%s/%s' % (str(num), str(num_primer))

    z_partial = data_primer[:max_ticks, :max_ticks]
    z_full = data[:max_ticks, :max_ticks]
    text = [
        [parse_labels(z_full[i, j], z_partial[i, j]) for j in range(z_full.shape[1])] for i in range(z_full.shape[0])
    ]

    z_partial_out: list[list[int]] = z_partial.tolist()
    z_full_out: list[list[int]] = z_full.tolist()
    return (z_partial_out, z_full_out, text, str1, str2)


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


def generate_motifb64(seq: str, row: dict) -> tuple:
    highlight = list(map(int, str(row['repetition_index']).split('_')))
    # print(f"{highlight=}") -> [1, 2]
    sequence, _ = highlight_subpart(seq, highlight)
    motif_name = row['motif_name']
    motif_name_part1 = f'{motif_name.replace("/", "_")}'
    motif_name_part2 = f'{",".join(map(str, highlight)) if highlight is not None else "mot"}'
    motif_name_long = f'{motif_name_part1}_{motif_name_part2}'
    motif_clean = re.sub(r'[^\w_]', '', motif_name_long)
    motif_id = motif_clean.rsplit('_', 1)[0]
    motif_clean_id = motif_id if highlight == [1] else motif_clean  # trick to solve static html
    # motif_clean_id sucks... and unfortunatelly it is used as module_id in json

    a1 = row['allele1']
    a2 = row['allele2']
    conf_total = float_to_str(row['confidence'], percents=True)
    conf_a1 = float_to_str(row['conf_allele1'], percents=True)
    conf_a2 = float_to_str(row['conf_allele2'], percents=True)
    if (a1 == 'B' and a2 == 'B') or (a1 == 0 and a2 == 0):
        result = f'BG {conf_total}'
    else:
        result = f'{str(a1):2s} ({conf_a1}) {str(a2):2s} ({conf_a2}) total {conf_total}'

    alignment = f"{motif_name}/alignments.html"

    return (motif_clean_id, highlight, motif_id, motif_name, result, alignment, sequence)


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


def generate_result_line(
    motif: Motif, predicted: tuple[str | int, str | int], confidence: tuple[float | str, ...],
    qual_num: int, primer_num: int, filt_num: int, module_number: int,
    qual_annot: list[Annotation] | None = None,
    second_module_number: int | None = None
) -> dict:
    """
    Generate result line from the template string.
    :param motif_class: Motif - motif class
    :param predicted: tuple[str, str] - predicted alleles (number or 'B'/'E')
    :param confidence: tuple[7x float/str] - confidences of prediction
    :param qual_num: int - number of reads with both primers
    :param primer_num: int - number of reads with exactly one primer
    :param filt_num: int - number of filtered out reads (no primers, many errors, ...)
    :param module_number: int - module number in motif
    :param qual_annot: list[Annotation] - list of quality annotations for error and number of reads
    :param second_module_number: int/None - second module number in motif
    :return: dict - result dictionary
    """
    # setup motif info
    start, end = motif.get_location_subpart(module_number)
    motif_seq = motif.module_str(module_number)
    repetition_index: int | str = module_number
    if second_module_number is not None:
        _, end = motif.get_location_subpart(second_module_number)
        motif_seq = ','.join([motif.module_str(i) for i in range(module_number, second_module_number + 1)])
        repetition_index = f'{module_number}_{second_module_number}'

    reads_a1: int | str
    reads_a2: int | str
    indels_rel: float | str
    indels_rel1: float | str
    indels_rel2: float | str
    mismatches_rel: float | str
    mismatches_rel1: float | str
    mismatches_rel2: float | str
    # get info about errors and number of reads from quality annotations if provided
    reads_a1 = reads_a2 = '---'
    indels_rel = mismatches_rel = '---'
    indels_rel1 = mismatches_rel1 = '---'
    indels_rel2 = mismatches_rel2 = '---'
    if qual_annot is not None:
        # get info about number of reads
        a1 = int(predicted[0]) if isinstance(predicted[0], int) else None
        a2 = int(predicted[1]) if isinstance(predicted[1], int) else None
        reads_a1 = 0 if a1 is None else len([a for a in qual_annot if a.module_repetitions[module_number] == a1])
        reads_a2 = 0 if a2 is None else len([a for a in qual_annot if a.module_repetitions[module_number] == a2])

        # get info about errors
        errors = [a.get_module_errors(motif, module_number) for a in qual_annot]
        errors_a1 = [a.get_module_errors(motif, module_number) for a in qual_annot
                     if a.module_repetitions[module_number] == a1]
        errors_a2 = [a.get_module_errors(motif, module_number) for a in qual_annot
                     if a.module_repetitions[module_number] == a2]
        assert len([l for i, _, l in errors if l == 0]) == 0

        # extract error metrics
        indels_rel, mismatches_rel = errors_per_read(errors, relative=True)
        indels_rel1, mismatches_rel1 = errors_per_read(errors_a1, relative=True)
        indels_rel2, mismatches_rel2 = errors_per_read(errors_a2, relative=True)

    return {
        'motif_name': motif.name, 'motif_nomenclature': motif.nomenclature, 'motif_sequence': motif_seq,
        'chromosome': motif.chrom, 'start': start, 'end': end,
        'allele1': predicted[0], 'allele2': predicted[1],
        'confidence': confidence[0], 'conf_allele1': confidence[1], 'conf_allele2': confidence[2],
        'reads_a1': reads_a1, 'reads_a2': reads_a2,
        'indels': indels_rel, 'indels_a1': indels_rel1, 'indels_a2': indels_rel2,
        'mismatches': mismatches_rel, 'mismatches_a1': mismatches_rel1, 'mismatches_a2': mismatches_rel2,
        'quality_reads': qual_num, 'one_primer_reads': primer_num, 'filtered_reads': filt_num,
        'conf_background': confidence[3] if len(confidence) > 3 else '---',
        'conf_background_all': confidence[4] if len(confidence) > 4 else '---',
        'conf_extended': confidence[5] if len(confidence) > 5 else '---',
        'conf_extended_all': confidence[6] if len(confidence) > 6 else '---',
        'repetition_index': repetition_index
    }


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


def generate_row(sequence: str, result: dict) -> tuple:
    """
    Generate rows of a summary table in html report.
    :param sequence: str - motif sequence
    :param result: pd.Series - result row to convert to table
    :param postfilter: PostFilter - postfilter dict from config
    :return: str - html string with rows of the summary table
    """
    highlight = list(map(int, str(result['repetition_index']).split('_')))
    sequence, _subpart = highlight_subpart(sequence, highlight)

    # shorten sequence:
    keep = 10
    first = sequence.find(',')
    last = sequence.rfind(',')
    smaller_seq = sequence if first == -1 else '...' + sequence[first - keep:last + keep + 1] + '...'

    # fill templates:
    updated_result = {
        'conf_allele1': float_to_str(result['conf_allele1'], percents=True),
        'conf_allele2': float_to_str(result['conf_allele2'], percents=True),
        'confidence': float_to_str(result['confidence'], percents=True),
        'motif_nomenclature': smaller_seq,
        'indels': float_to_str(result['indels'], decimals=2),
        'mismatches': float_to_str(result['mismatches'], decimals=2),
        'indels_a1': float_to_str(result['indels_a1'], decimals=2),
        'mismatches_a1': float_to_str(result['mismatches_a1'], decimals=2),
        'indels_a2': float_to_str(result['indels_a2'], decimals=2),
        'mismatches_a2': float_to_str(result['mismatches_a2'], decimals=2)
    }
    # return ROW_STRING.format(**{**result, **updated_result})
    motif_name = result['motif_name']
    motif_nomenclature = updated_result['motif_nomenclature']
    allele1 = result['allele1']
    conf_allele1 = updated_result['conf_allele1']
    reads_a1 = result['reads_a1']
    indels_a1 = updated_result['indels_a1']
    mismatches_a1 = updated_result['mismatches_a1']

    allele2 = result['allele2']
    conf_allele2 = updated_result['conf_allele2']
    reads_a2 = result['reads_a2']
    indels_a2 = updated_result['indels_a2']
    mismatches_a2 = updated_result['mismatches_a2']

    confidence = updated_result['confidence']
    indels = updated_result['indels']
    mismatches = updated_result['mismatches']
    quality_reads = result['quality_reads']
    one_primer_reads = result['one_primer_reads']

    row_tuple = (
        motif_name, motif_nomenclature,
        allele1, conf_allele1, reads_a1, indels_a1, mismatches_a1,
        allele2, conf_allele2, reads_a2, indels_a2, mismatches_a2,
        confidence, indels, mismatches, quality_reads, one_primer_reads
    )

    return row_tuple
    # return (result, updated_result)


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

    def modules_str(self, include_flanks: bool = False) -> str:
        """
        Returns string representation of modules
        :param include_flanks: bool - include flank modules?
        :return: String representation of modules
        """
        if include_flanks:
            return ''.join([f'{seq}[{num}]' for seq, num in self.modules])
        return ''.join([f'{seq}[{num}]' for seq, num in self.modules[1:-1]])

    def module_str(self, module_number: int) -> str:
        """
        Returns string representation of modules
        :param module_number: int - module number
        :return: String representation of modules
        """
        seq, num = self.modules[module_number]
        return f'{seq}[{num}]'

    def get_location_subpart(self, index: int) -> tuple[int, int]:
        """
        Returns the chromosome location of a subpart of a motive
        :param index: int - index of a module
        :return: start and end location of the subpart
        """
        start = self.start
        for module in self.modules[1: index]:
            seq, rep = module
            start += len(seq) * rep

        return start, start + len(self.modules[index][0]) * self.modules[index][1]


class Annotation:
    def __init__(self, row: pd.Series):
        self.states = row["modules"]
        self.mismatches_string = row["mismatches_str"]
        self.module_repetitions = list(map(int, row["module_repetitions"].split(",")))
        self.module_nomenclatures = list(row["module_nomenclatures"].split(","))

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

    def is_annotated_right(self) -> bool:
        """
        Is it annotated in a way that it is interesting?
        More than one module annotated + modules are not missing in the middle.
        :return: bool - annotated right?
        """

        # remove those that starts/ends in background but don't have a neighbour module
        starts_background = self.states[0] in '_-'
        ends_background = self.states[-1] in '_-'
        if starts_background and self.module_repetitions[0] == 0:
            return False
        if ends_background and self.module_repetitions[-1] == 0:
            return False

        # remove those with jumping modules
        started = False
        ended = False
        for repetition in self.module_repetitions:
            if repetition > 0:
                started = True
                if ended:
                    return False
            if repetition == 0 and started:
                ended = True

        # pass?
        return True

    def get_str_repetitions(self, index_str: int) -> tuple[bool, int] | None:
        """
        Get the number of str repetitions for a particular index.
        :param index_str: int - index of a str
        :return: (bool, int) - closed?, number of str repetitions
        """
        if self.is_annotated_right():
            primer1 = index_str > 0 and self.module_repetitions[index_str - 1] > 0
            primer2 = index_str + 1 < len(self.module_repetitions) and self.module_repetitions[index_str + 1] > 0
            if primer1 or primer2:
                return primer1 and primer2, self.module_repetitions[index_str]
        return None
