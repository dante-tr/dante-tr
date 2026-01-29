from __future__ import annotations  # mute typechecking of classes declared later than used
from jinja2 import Environment, FileSystemLoader  # type: ignore
import pandas as pd
import numpy as np
import shutil
import enum
import json
import sys
import re
import os

MOTIF_COLUMN_NAME = 'motif'
MIN_FLANK_LEN = 3
MIN_REP_LEN = 3
MIN_REP_CNT = 1
MAX_ABS_ERROR = None
MAX_REL_ERROR = 1.0
MAX_REPETITIONS = 40
BASE_MAPPING = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'M': '[AC]', 'R': '[AG]', 'W': '[AT]', 'S': '[CG]', 'Y': '[CT]', 'K': '[GT]',
    'V': '[ACG]', 'H': '[ACT]', 'D': '[AGT]', 'B': '[CGT]',
    'N': '[ACGT]'
}


def main() -> None:
    input_tsvs = ["457-2025_WGS_FAME3/annotations.tsv"]
    input_jsons = ["457-2025_WGS_FAME3/data.json"]
    output_dir = "."
    output_files = [f"{output_dir}/alignments/FAME3.html"]
    is_male = True

    for (tsv_file, json_file, output_file) in zip(input_tsvs, input_jsons, output_files):
        # os.makedirs(motif_dir, exist_ok=True)  # create output directory for alignments?
        write_alignment_html(json_file, tsv_file, output_file, is_male)
    copy_includes(output_dir)


def write_alignment_html(json_file, tsv_file, output_file, is_male) -> None:
    with open(json_file, "r") as f:
        data_json = json.load(f)

    df = pd.read_csv(tsv_file, sep='\t')

    motif = create_motif(df, is_male)
    mt = motif.dir_name()
    seq = motif.modules_str(include_flanks=True)
    motif_desc = motif.name
    data = []
    for motif_data in data_json["motifs"]:
        assert motif_data["motif_id"] == motif_desc, "Some incompatible input"

        for module in motif_data["modules"]:
            mod_id = module["id"][0]
            a1 = int(module["allele_1"][0]) if isinstance(module["allele_1"][0], int) else None
            a2 = int(module["allele_2"][0]) if isinstance(module["allele_2"][0], int) else None
            fastas = gen_single_fastas(tsv_file, is_male, mod_id, a1, a2)
            data.append(format_single_msas(mod_id, seq, motif_desc, a1, a2, fastas))

        for phasing in motif_data["phasings"]:
            md1_id = phasing["ids"][0]
            md1_id = phasing["ids"][1]
            fastas = gen_phased_fastas(tsv_file, is_male, md1_id, md1_id)
            data.append(format_phased_msas(md1_id, md1_id, seq, motif_desc, fastas))

    script_dir = os.path.dirname(sys.argv[0]) + "/templates"
    env = Environment(loader=FileSystemLoader([script_dir]))
    template = env.get_template("alignments_template.html")
    output = template.render(sample=mt, motif_desc=motif_desc, data2=data)
    with open(output_file, "w") as f:
        f.write(output)

    return


def gen_single_fastas(input_tsv, is_male, m_idx, a1, a2) -> list[str]:
    result: list[str] = [
        gen_msa(input_tsv, "spanning", is_male, m_idx, a1, None, None, False),
        gen_msa(input_tsv, "spanning", is_male, m_idx, a2, None, None, False),
        gen_msa(input_tsv, "spanning", is_male, m_idx, None, None, None, False),
        gen_msa(input_tsv, "flanking", is_male, m_idx, None, None, None, False),
        gen_msa(input_tsv, "flanking_left", is_male, m_idx, None, None, None, False),
        gen_msa(input_tsv, "flanking_right", is_male, m_idx, None, None, None, True),
    ]
    return result


def gen_phased_fastas(input_tsv, is_male, m1_idx, m2_idx) -> list[str]:
    result: list[str] = [
        gen_msa(input_tsv, "two_good", is_male, m1_idx, None, m2_idx, None, False),
        gen_msa(input_tsv, "one_good", is_male, m1_idx, None, m2_idx, None, False),
        gen_msa(input_tsv, "one_good_left", is_male, m1_idx, None, m2_idx, None, False),
        gen_msa(input_tsv, "one_good_right", is_male, m1_idx, None, m2_idx, None, True),
    ]
    return result


def format_single_msas(m_, seq, motif_desc, a1, a2, fastas) -> tuple:
    data = []
    highlight = list(map(int, str(m_).split('_')))
    sequence, _subpart = highlight_subpart(seq, highlight)
    motif_name_part1 = f'{motif_desc.replace("/", "_")}'
    motif_name_part2 = f'{",".join(map(str, highlight)) if highlight is not None else "mot"}'
    motif_name = f'{motif_name_part1}_{motif_name_part2}'
    motif_clean = re.sub(r'[^\w_]', '', motif_name)
    motif_id = motif_clean.split('_')[0]

    if a1 is not None and a1 > 0:
        name = f'{motif_clean}_{str(a1)}'
        display_text = f'Allele 1 ({str(a1):2s}) alignment visualization'
        fasta2 = fastas[0]
        seq_logo = 'true'
        data.append((name, display_text, fasta2, seq_logo))

    if a2 is not None and a2 != a1 and a2 != 0:
        name = f'{motif_clean}_{str(a2)}'
        display_text = f'Allele 2 ({str(a2):2s}) alignment visualization'
        fasta2 = fastas[1]
        seq_logo = 'true'
        data.append((name, display_text, fasta2, seq_logo))

    name = motif_clean
    display_text = 'Spanning reads alignment visualization'
    fasta2 = fastas[2]
    seq_logo = 'true'
    data.append((name, display_text, fasta2, seq_logo))

    name = motif_clean + '_filtered'
    display_text = 'Partial reads alignment visualization'
    fasta2 = fastas[3]
    seq_logo = 'false'
    data.append((name, display_text, fasta2, seq_logo))

    name = motif_clean + '_filtered_left'
    display_text = 'Left flank reads alignment visualization'
    fasta2 = fastas[4]
    seq_logo = 'true'
    data.append((name, display_text, fasta2, seq_logo))

    name = motif_clean + '_filtered_right'
    display_text = 'Right flank reads alignment visualization'
    fasta2 = fastas[5]
    seq_logo = 'true'
    data.append((name, display_text, fasta2, seq_logo))

    return (sequence, motif_id, data)


def format_phased_msas(m1, m2, seq, motif_desc, fastas2) -> tuple:
    data = []

    suffix = f'{m1}_{m2}'
    highlight = list(map(int, str(suffix).split('_')))
    sequence, _subpart = highlight_subpart(seq, highlight)
    motif_name_part1 = f'{motif_desc.replace("/", "_")}'
    motif_name_part2 = f'{",".join(map(str, highlight)) if highlight is not None else "mot"}'
    motif_name = f'{motif_name_part1}_{motif_name_part2}'
    motif_clean = re.sub(r'[^\w_]', '', motif_name)
    motif_id = motif_clean.split('_')[0]

    name = motif_clean
    display_text = 'Spanning reads alignment visualization'
    fasta2 = fastas2[0]
    seq_logo = 'true'
    data.append((name, display_text, fasta2, seq_logo))

    name = motif_clean + '_filtered'
    display_text = 'Partial reads alignment visualization'
    fasta2 = fastas2[1]
    seq_logo = 'true'
    data.append((name, display_text, fasta2, seq_logo))

    name = motif_clean + '_filtered_left'
    display_text = 'Left flank reads alignment visualization'
    fasta2 = fastas2[2]
    seq_logo = 'true'
    data.append((name, display_text, fasta2, seq_logo))

    name = motif_clean + '_filtered_right'
    display_text = 'Right flank reads alignment visualization'
    fasta2 = fastas2[3]
    seq_logo = 'true'
    data.append((name, display_text, fasta2, seq_logo))

    return (sequence, motif_id, data)


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

    def __init__(self, motif: str, name: str | None = None) -> None:
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
        self.monoallelic: bool = False

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

    def module_str(self, module_number: int) -> str:
        """
        Returns string representation of modules
        :param module_number: int - module number
        :return: String representation of modules
        """
        seq, num = self.modules[module_number]
        return f'{seq}[{num}]'

    def dir_name(self) -> str:
        """
        Returns possible directory name of the motif.
        :return: str - directory name for the motif
        """
        return self.name

    def get_repeating_modules(self) -> list[tuple[int, str, int]]:
        """
        Returns list of modules with more than one repetition.
        :return: List of tuples containing index, sequence, and repetition count.
        """
        return [(int(i), str(seq), int(num)) for i, (seq, num) in enumerate(self.modules) if num > 1]

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


def gen_msa(
    input_tsv: str,
    filter_type: str, is_male: bool,
    m1: int,
    a1: int | None,
    m2: int | None,
    a2: int | None,
    right_align: bool
) -> str:
    print(input_tsv, filter_type, is_male, m1, a1, m2, a2, right_align)
    df = pd.read_csv(input_tsv, sep='\t')

    annotations = select_relevant_annotations(df, is_male, m1, m2, filter_type, a1, a2)
    print(len(annotations))

    # apply cutoff
    motif = create_motif(df, is_male)
    annotations = [a.get_shortened_annotation(5, motif) for a in annotations]

    # setup alignments
    alignments, states = setup_alignments(annotations)

    # sort according to motif count:
    left_flank = np.array([-(ann.module_bases[0] + ann.left_flank_len) for ann in annotations])
    left_flank_exist = np.array([-(ann.module_repetitions[0]) if not right_align else 0 for ann in annotations])
    if m2 is not None:
        # sorting first with 1st allele then with second
        reps1 = np.array([-ann.module_bases[m1] for ann in annotations])
        reps2 = np.array([-ann.module_bases[m2] for ann in annotations])
        # sort by existence of left flank, first allele, second, left flank len.
        sort_inds = np.lexsort((left_flank, reps2, reps1, left_flank_exist))
    else:
        reps = np.array([-ann.module_bases[m1] for ann in annotations])
        sort_inds = np.lexsort((left_flank, reps, left_flank_exist))
    annotations = np.array(annotations)[sort_inds]
    alignments = list(np.array(alignments)[sort_inds])

    # for every alignment, shift the left flank right
    end = get_left_flank(states)
    if end != -1:
        for i in range(len(alignments)):
            alignments[i] = move_right(alignments[i], 0, end)

    # for every alignment, shift the first module right
    start0, end0 = get_range(states, '0')
    if start0 != -1:
        for i in range(len(alignments)):
            alignments[i] = move_right(alignments[i], start0, end0)

    # in addition, those that have only '_' in state '0' (missing left flank), shift right also '1' state
    start1, end1 = get_range(states, '1')
    first_zero_idx = len(alignments)
    for i in range(len(alignments)):
        if start1 != -1 and (start0 == -1 or alignments[i][start0:end0].count('_') == end0 - start0 or right_align):
            first_zero_idx = min(first_zero_idx, i)
            alignments[i] = move_right(alignments[i], start1, end1)

    # add empty line if we have some alignments without left flank
    annot_names = [annot.read_id for annot in annotations]
    if first_zero_idx != len(alignments) and not right_align:
        alignments = alignments[:first_zero_idx] + ['_' * len(alignments[0])] + alignments[first_zero_idx:]
        annot_names = annot_names[:first_zero_idx] + ['empty_line'] + annot_names[first_zero_idx:]

    string_tmp = []
    for annot_name, align in zip(annot_names, alignments):
        string_tmp.append(f">{annot_name}\n{align}\n")
    string = "".join(string_tmp)

    return string


def select_relevant_annotations(df, is_male, m1, m2, filter_type, a1, a2):
    # df["cls_left"] = df["module_classes"].apply(lambda x: x.split(",")[0])
    # df["cls"] = df["module_classes"].apply(lambda x: x.split(",")[1])
    # df["cls_right"] = df["module_classes"].apply(lambda x: x.split(",")[2])

    motif = create_motif(df, is_male)
    annotations = create_annotations(df, motif)

    postfilter = PostFilter()
    anns_spanning, rest = postfilter.get_filtered(motif, annotations, m1, both_primers=True)
    anns_flanking, anns_filtered = postfilter.get_filtered(motif, rest, m1, both_primers=False)
    del rest

    anns_flank_lt = [a for a in anns_flanking if a.module_bases[0] > 0]
    anns_flank_rt = [a for a in anns_flanking if a.module_bases[-1] > 0]

    if m2 is not None:
        mod_nums = [m1, m2]
        anns_good_two, _filtered = postfilter.get_filtered_list(motif, annotations, mod_nums, both_primers=[True, True])
        _left_good, _left_bad = postfilter.get_filtered_list(motif, _filtered, mod_nums, both_primers=[False, True])
        _right_good, anns_0good = postfilter.get_filtered_list(motif, _left_bad, mod_nums, both_primers=[True, False])
        anns_good_one = _left_good + _right_good
        del _filtered, _left_good, _left_bad, _right_good

        anns_1good_lt = [a for a in anns_good_one if a.module_bases[0] > 0]
        anns_1good_rt = [a for a in anns_good_one if a.module_bases[-1] > 0]

    if filter_type == "spanning":
        annotations = select_annotation(anns_spanning, a1, a2, m1, m2)
        # annotations = anns_spanning
    elif filter_type == "flanking":
        annotations = select_annotation(anns_flanking, a1, a2, m1, m2)
        # annotations = anns_flanking
    elif filter_type == "flanking_left":
        annotations = select_annotation(anns_flank_lt, a1, a2, m1, m2)
        # annotations = anns_flank_lt
    elif filter_type == "flanking_right":
        annotations = select_annotation(anns_flank_rt, a1, a2, m1, m2)
        # annotations = anns_flank_rt
    elif filter_type == "two_good":
        annotations = select_annotation(anns_good_two, a1, a2, m1, m2)
        # annotations = anns_good_two
    elif filter_type == "one_good":
        annotations = select_annotation(anns_good_one, a1, a2, m1, m2)
        # annotations = anns_good_one
    elif filter_type == "one_good_left":
        annotations = select_annotation(anns_1good_lt, a1, a2, m1, m2)
        # annotations = anns_1good_lt
    elif filter_type == "one_good_right":
        annotations = select_annotation(anns_1good_rt, a1, a2, m1, m2)
        # annotations = anns_1good_rt

    print(len(annotations))
    return annotations


class ChromEnum(enum.Enum):
    X = 'X'
    Y = 'Y'
    NORM = 'NORM'


def chrom_from_string(chrom_str: str) -> ChromEnum:
    """
    Converts a string to a ChromEnum object.
    :param chrom_str: str - the string to convert to a ChromEnum object
    :return ChromEnum - enum object representing the chromosome
    """
    return (
        ChromEnum.X if chrom_str in ['chrX', 'NC_000023'] else
        ChromEnum.Y if chrom_str in ['chrY', 'NC_000024'] else
        ChromEnum.NORM
    )


def select_annotation(annotations, allele1, allele2, index_rep1, index_rep2):
    print(len(annotations))
    if allele1 is None:
        return annotations
    if allele2 is None or index_rep2 is None:
        return [a for a in annotations if a.module_repetitions[index_rep1] == allele1]

    result = []
    for a in annotations:
        if a.module_repetitions[index_rep1] == allele1 and a.module_repetitions[index_rep2] == allele2:
            result.append(a)
    return result


def get_left_flank(states) -> int:
    """
    Get range of a left flank.
    :return: int - (one after) nd of the left flank before module '0'
    """
    for i, state in enumerate(states):
        if state != '-':
            return i
    return -1


def move_right(alignment: str, start: int, end: int) -> str:
    """
    Shift first part of the alignment to the right.
    :param alignment: str - alignment of the read
    :param start: int - start idx for shift
    :param end: int - one after end idx for a shift
    :return: str - alignment, where first part is shifted to right
    """
    align_part = alignment[start:end]
    # find last empty:
    idx = 0
    for idx in reversed(range(len(align_part))):
        if align_part[idx] != '_':
            break
    idx += 1

    # return shifted alignment
    return alignment[:start] + ('_' * (len(align_part) - idx)) + align_part[:idx] + alignment[end:]


def get_range(states, symbol: str) -> tuple[int, int]:
    """
    Get range of a module
    :param symbol: str - symbol for state to get range for
    :return: int, int - start and (one after) end range coordinates
    """
    try:
        first_idx = states.index(symbol)
        last_idx = len(states) - states[-1::-1].index(symbol) - 1
        return first_idx, last_idx + 1
    except ValueError:
        return -1, -1


class Annotation:
    """
    Encapsulate sequence of states from HMM and provide its readable representation and filters
    """

    def __init__(
        self, read_id: str, mate_order: int, read_seq: str, expected_seq: str,
        states: str, probability: float, motif: Motif
    ):
        """
        :param read_id: str - read ID
        :param read_seq: str - read sequence
        :param mate_order: int - mate order (0 - unpaired, 1 - left pair, 2 - right pair)
        :param expected_seq: str - expected sequence as in motif
        :param states: str - sequence of states (numbers of modules)
        :param probability: Probability of generating sequence by the most likely sequence of HMM states
        :param motif: Sequence of tuples (sequence, repeats) as specified by user
        """

        # Store arguments into instance variables
        self.read_id = read_id
        self.mate_order = mate_order
        self.read_seq = read_seq
        self.expected_seq = expected_seq
        self.states = states
        self.probability = probability
        self.n_modules = len(motif.modules)

        # Calculate insertion/deletion/mismatch string
        self.mismatches_string = self.__get_errors()

        # Calculate number of insertions, deletions and normal bases
        self.n_insertions = self.mismatches_string.count('I')
        self.n_deletions = self.mismatches_string.count('D')
        self.n_mismatches = self.mismatches_string.count('M')

        # Number of STR motif repetitions and sequences of modules
        self.module_bases = self.__get_bases_per_module()
        self.module_repetitions = self.__get_module_repetitions(motif)
        self.module_sequences = self.__get_module_sequences()

        # get left flank length
        self.left_flank_len = self.__get_left_flank()

    def __str__(self) -> str:
        """
        Return the annotation.
        :return: str - annotation
        """
        return '\n'.join([f'{self.read_id} {str(self.module_bases)} {str(self.module_repetitions)}', self.read_seq,
                          self.expected_seq, self.states, self.mismatches_string])

    def __get_errors(self) -> str:
        """
        Count errors in annotation and the error line.
        :return: str - error line
        """
        err_line = []
        for exp, read in zip(self.expected_seq.upper(), self.read_seq.upper()):
            if exp == '-' or read in BASE_MAPPING.get(exp, ''):
                err_line.append('_')
            elif read == '_':
                err_line.append('D')
            elif exp == '_':
                err_line.append('I')
            else:
                err_line.append('M')

        return ''.join(err_line)

    def __get_bases_per_module(self) -> tuple[int, ...]:
        """
        List of integers, each value corresponds to number of bases of input sequence that were generated by the module
        :return: Number of bases generated by each module
        """
        # Count the module states
        return tuple(self.states.count(chr(ord('0') + i)) for i in range(self.n_modules))

    def __get_left_flank(self) -> int:
        """
        Get length of a left flank.
        :return: int - number of bases of left flank before module '0' (usually module '0' is still left flank)
        """
        for i, state in enumerate(self.states):
            if state != '-':
                return i
        return len(self.states)

    def __get_module_repetitions(self, motif: Motif) -> tuple[int, ...]:
        """
        List of integers, each value corresponds to number of repetitions of module in annotation
        :return: Number of repetitions generated by each module
        """
        # Count the module states
        repetitions = self.__get_bases_per_module()

        # Divide by the module length where applicable
        # TODO: this is not right for grey ones, where only closed ones should be counted, so round is not right.
        return tuple(
            1 if reps == 1 and cnt > 0 else round(cnt / len(seq))
            for (seq, reps), cnt in zip(motif.modules, repetitions)
        )

    def __get_module_sequences(self) -> tuple[str, ...]:
        """
        List of sequences, each per module
        :return: list(str)
        """
        sequences = [''] * self.n_modules
        for i in range(self.n_modules):
            state_char = chr(ord('0') + i)
            first = self.states.find(state_char)
            if first > -1:
                last = self.states.rfind(state_char)
                sequences[i] = self.read_seq[first:last + 1]
        return tuple(sequences)

    def get_module_errors(self, motif: Motif, module_num: int, overhang: int | None = None) -> tuple[int, int, int]:
        """
        Get the number of insertions and deletions or mismatches in a certain module.
        If overhang is specified, look at specified number of bases around the module as well.
        :param module_num: int - 0-based module number to count errors
        :param overhang: int - how long to look beyond module, if None, one length of STR module
        :return: int, int, int - number of insertions and deletions, mismatches, length of the interval
        """
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

    def has_less_errors(self, max_errors: float | int | None, relative=False) -> bool:
        """
        Check if this annotation has fewer errors than max_errors.
        Make it relative to the annotated length if relative is set.
        :param max_errors: int/float - number of max_errors (relative if relative is set)
        :param relative: bool - if the errors are relative to the annotated length
        :return: bool - True if the number of errors is less than allowed
        """
        errors = self.n_deletions + self.n_insertions + self.n_mismatches

        if max_errors is None or errors == 0:
            return True

        if relative:
            return errors / float(sum(self.module_bases)) <= max_errors
        return errors <= max_errors

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

    @staticmethod
    def find_with_regex(read_sequence: str, motif_sequence: str, search_pos: int = 0) -> int:
        """
        Find the first occurrence of a motif sequence in the read sequence using regular expressions.
        :param read_sequence: The sequence to search in.
        :param motif_sequence: The motif sequence (as a regex) to search for.
        :param search_pos: The position to start the search from.
        :return: int - The start position of the first occurrence of the motif sequence. Returns -1 if not found.
        """
        # convert motif sequence to regex
        motif_regex = ''.join(BASE_MAPPING[char] for char in motif_sequence)

        # compile the regular expression pattern
        pattern = re.compile(motif_regex)

        # search for the pattern in the read sequence starting from search_pos
        match = pattern.search(read_sequence, search_pos)

        # return the start position if a match is found, else return -1
        return match.start() if match else -1

    def get_nomenclature(
        self, motif: Motif, index_rep: int | None = None, index_rep2: int | None = None, include_flanking: bool = True
    ) -> str:
        """
        Get HGVS nomenclature.
        :param index_rep: int - index of the first repetition (None if include all)
        :param index_rep2: int - index of the second repetition (None if include all)
        :param include_flanking: boolean - include flanking regions (i.e. first and last module)
        :return: str - HGVS nomenclature string
        """
        # prepare data
        if index_rep is not None:
            if index_rep2 is not None:
                data = zip(
                    [self.module_repetitions[index_rep], self.module_repetitions[index_rep2]],
                    [motif[index_rep], motif[index_rep2]],
                    [self.module_sequences[index_rep], self.module_sequences[index_rep2]]
                )
            else:
                data = zip(
                    [self.module_repetitions[index_rep]],
                    [motif[index_rep]],
                    [self.module_sequences[index_rep]]
                )
        elif include_flanking:
            data = zip(self.module_repetitions, motif.modules, self.module_sequences)
        else:
            data = zip(self.module_repetitions[1:-1], motif.modules[1:-1], self.module_sequences[1:-1])

        # iterate and build the nomenclature string
        nomenclatures = []
        for repetitions, (motif_sequence, _), read_sequence in data:
            nomenclature = self.build_nomenclature_string(repetitions, motif_sequence, read_sequence)
            nomenclatures.append(nomenclature)

        return '\t'.join(nomenclatures)

    def build_nomenclature_string(self, repetitions, motif_sequence, read_sequence) -> str:
        nomenclature = ''
        if repetitions == 1:
            if len(read_sequence) > 0:
                nomenclature += f'{read_sequence}[1]'
            return nomenclature

        reps = 0
        search_pos = 0
        found_rep_seq = ''
        while True:
            search_found = self.find_with_regex(read_sequence, motif_sequence, search_pos)
            if search_found == search_pos:
                # setup current rep. sequence
                if reps == 0:
                    found_rep_seq = read_sequence[search_found:search_found + len(motif_sequence)]

                if read_sequence[search_found:search_found + len(motif_sequence)] == found_rep_seq:
                    # regular continuation
                    reps += 1
                else:
                    # interruption, but in line with searched motif
                    nomenclature += f'{found_rep_seq}[{reps}]'
                    found_rep_seq = read_sequence[search_found:search_found + len(motif_sequence)]
                    reps = 1
            elif search_found == -1:
                # the end, we did not find any other STRs
                if reps > 0:
                    nomenclature += f'{found_rep_seq}[{reps}]'
                if len(read_sequence[search_pos:]) > 0:
                    nomenclature += f'{read_sequence[search_pos:]}[1]'
                break
            else:
                # some interruption
                if reps > 0:
                    nomenclature += f'{found_rep_seq}[{reps}]'
                if len(read_sequence[search_pos:search_found]) > 0:
                    nomenclature += f'{read_sequence[search_pos:search_found]}[1]'
                found_rep_seq = read_sequence[search_found:search_found + len(motif_sequence)]
                reps = 1
            # update search pos and iterate
            search_pos = search_found + len(motif_sequence)
        return nomenclature

    def get_shortened_annotation(self, shorten_length: int, motif: Motif) -> Annotation:
        # search for start
        start = -1
        for i in range(len(self.states)):
            if self.states[i] != '-':
                start = i
                break
        start = max(start - shorten_length, 0)

        # search for end
        end = -1
        for i in range(len(self.states) - 1, -1, -1):
            if self.states[i] != '-':
                end = i
                break
        end = min(end + 1 + shorten_length, len(self.states))  # +1 for use as list range

        # return shortened Annotation
        return Annotation(
            self.read_id, self.mate_order, self.read_seq[start:end], self.expected_seq[start:end],
            self.states[start:end], self.probability, motif
        )


def create_motif(df: pd.DataFrame, is_male: bool) -> Motif:
    motif_str = df[MOTIF_COLUMN_NAME].iloc[0]
    name = None if 'name' not in df.columns or df.iloc[0]['name'] in ['None', ''] else df.iloc[0]['name']
    motif_class = Motif(motif_str, name)
    motif_class.monoallelic = is_male and chrom_from_string(motif_class.chrom) in [ChromEnum.X, ChromEnum.Y]
    return motif_class


def create_annotations(df: pd.DataFrame, motif: Motif) -> list[Annotation]:
    annotations: list[Annotation] = []
    for _, row in df.iterrows():
        ann = Annotation(
            row['read_id'], row['mate_order'], row['read'], row['reference'],
            row['modules'], row['log_likelihood'], motif
        )
        annotations.append(ann)
    return annotations


class PostFilter:
    """
    Class that encapsulates post-filtering.
    """

    def __init__(self):
        self.min_flank_len = MIN_FLANK_LEN
        self.min_rep_len = MIN_REP_LEN
        self.min_rep_cnt = MIN_REP_CNT
        self.max_rel_error = MAX_REL_ERROR
        self.max_abs_error = MAX_ABS_ERROR

    def quality_annotation(self, motif: Motif, ann: Annotation, module_number: int, both_primers: bool = True) -> bool:
        """
        Is this annotation good?
        :param ann: Annotation - annotation to be evaluated
        :param module_number: int - module number
        :param both_primers: bool - do we require both primers to be present
        :return: bool - quality annotation?
        """
        is_right = ann.is_annotated_right()

        primers = ann.primers(module_number)
        has_primers = primers == 2 if both_primers else primers >= 1

        has_less_errors = (
            ann.has_less_errors(self.max_rel_error, relative=True)
            and ann.has_less_errors(self.max_abs_error, relative=False)
        )

        left_flank = sum(ann.module_bases[module_number + 1:]) >= self.min_flank_len
        right_flank = sum(ann.module_bases[:module_number]) >= self.min_flank_len
        has_flanks = left_flank and right_flank if both_primers else left_flank or right_flank

        has_repetitions = (
            ann.module_bases[module_number] >= self.min_rep_len
            and ann.module_repetitions[module_number] >= self.min_rep_cnt
        )

        _seq, reps = motif.modules[module_number]

        return is_right and has_primers and has_less_errors and has_flanks and (has_repetitions or reps == 1)

    def get_filtered_list(
        self, motif: Motif, annotations: list[Annotation],
        module_number: list[int], both_primers: list[bool] | None = None
    ) -> tuple[list[Annotation], list[Annotation]]:
        """
        Get filtered annotations (list of modules).
        :param annotations: list(Annotation) - annotations
        :param module_number: list(int) - module numbers
        :param both_primers: list(bool) or None - do we require both primers to be present
        :return: list(Annotation), list(Annotation) - quality annotations, non-quality annotations
        """
        # adjust input if needed
        if both_primers is None:
            both_primers = [True] * len(module_number)
        assert len(both_primers) == len(module_number)

        # filter annotations
        quality_annotations = [
            an for an in annotations
            if all((
                self.quality_annotation(motif, an, mn, both_primers=bp) for mn, bp in zip(module_number, both_primers)
            ))
        ]
        filtered_annotations = [an for an in annotations if an not in quality_annotations]

        return quality_annotations, filtered_annotations

    # TODO: is this just a specialized version of the previous method? Do we need this?
    def get_filtered(
        self, motif: Motif, annotations: list[Annotation], module_number: int, both_primers: bool = True
    ) -> tuple[list[Annotation], list[Annotation]]:
        """
        Get filtered annotations.
        :param annotations: list(Annotation) - annotations
        :param module_number: int - module number
        :param both_primers: bool - do we require both primers to be present
        :return: list(Annotation), list(Annotation) - quality annotations, non-quality annotations
        """
        # pick quality annotations
        quality_annotations = []
        filtered_annotations = []
        for an in annotations:
            if self.quality_annotation(motif, an, module_number, both_primers):
                quality_annotations.append(an)
            else:
                filtered_annotations.append(an)

        return quality_annotations, filtered_annotations


def setup_alignments(annotations: list[Annotation]) -> tuple[list[str], list[str]]:
    alignments = [''] * len(annotations)  # alignment strings
    align_inds = np.zeros(len(annotations), dtype=int)  # indices of annotations that were processed
    states = []  # has numbers of states in the final multiple alignment

    while True:
        # get minimal state:
        min_comp = (True, 'Z', -1)
        total_done = 0
        for i, (annot, ai) in enumerate(zip(annotations, align_inds)):
            if ai >= len(annot.states):
                total_done += 1
                continue
            state = annot.states[ai]
            comparator = (state != 'I', state, i)
            min_comp = min(comparator, min_comp)

        # if we have done every state, end:
        if total_done >= len(alignments):
            break

        states.append(min_comp[1])

        # now print all states, that are minimal:
        for i, (annot, ai) in enumerate(zip(annotations, align_inds)):
            if ai >= len(annot.states):
                alignments[i] += '_'  # put ends of each alignment to be of same length
                continue
            if annot.states[ai] == min_comp[1]:
                alignments[i] += annot.read_seq[ai]
                align_inds[i] += 1
            else:
                alignments[i] += '_'

    return alignments, states


def copy_includes(output_dir: str) -> None:
    include_dir = os.path.dirname(sys.argv[0]) + "/includes"
    os.makedirs(f'{output_dir}/includes', exist_ok=True)
    shutil.copy2(f'{include_dir}/msa.min.gz.js',            f'{output_dir}/includes/msa.min.gz.js')
    shutil.copy2(f'{include_dir}/plotly-2.14.0.min.js',     f'{output_dir}/includes/plotly-2.14.0.min.js')
    shutil.copy2(f'{include_dir}/jquery-3.6.1.min.js',      f'{output_dir}/includes/jquery-3.6.1.min.js')
    shutil.copy2(f'{include_dir}/datatables.min.js',        f'{output_dir}/includes/datatables.min.js')
    shutil.copy2(f'{include_dir}/styles.css',               f'{output_dir}/includes/styles.css')
    shutil.copy2(f'{include_dir}/w3.css',                   f'{output_dir}/includes/w3.css')
    shutil.copy2(f'{include_dir}/jquery.dataTables.css',    f'{output_dir}/includes/jquery.dataTables.css')


if __name__ == '__main__':
    main()
