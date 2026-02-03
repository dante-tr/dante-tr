from __future__ import annotations  # mute typechecking of classes declared later than used
from jinja2 import Environment, FileSystemLoader  # type: ignore
from typing import TypeAlias  # Any
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
MSA: TypeAlias = list[tuple[str, str]]


def main() -> None:
    input_tsvs = ["457-2025_WGS_FAME3/annotations.tsv", "432-2025_WGS_DM2/annotations.tsv"]
    input_jsons = ["457-2025_WGS_FAME3/data.json", "432-2025_WGS_DM2/data.json"]
    output_dir = "."
    output_files = [f"{output_dir}/alignments/FAME3.html", f"{output_dir}/alignments/DM2.html"]
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
            fastas = gen_single_fastas(tsv_file, mod_id, a1, a2)
            data.append(format_single_msas(mod_id, seq, motif_desc, a1, a2, fastas))

        for phasing in motif_data["phasings"]:
            md1_id = phasing["ids"][0]
            md2_id = phasing["ids"][1]
            fastas = gen_phased_fastas(tsv_file, md1_id, md2_id)
            data.append(format_phased_msas(md1_id, md1_id, seq, motif_desc, fastas))

    script_dir = os.path.dirname(sys.argv[0]) + "/templates"
    env = Environment(loader=FileSystemLoader([script_dir]))
    template = env.get_template("alignments_template.html")
    output = template.render(sample=mt, motif_desc=motif_desc, data2=data)
    with open(output_file, "w") as f:
        f.write(output)

    return


def gen_single_fastas(input_tsv, m_idx, all1, all2) -> list[str]:
    result: list[str] = [
        # These functions look similar, but are subtly different.
        # I think it is better to keep them separate, to make reasoning within them easier - less cases to keep in mind.
        # Also, their interface is simple - given pandas DataFrame in tsv and module, they return str representing MSA.
        gen_spann_al_msa(input_tsv, m_idx, all1),
        gen_spann_al_msa(input_tsv, m_idx, all2),
        gen_spanning_msa(input_tsv, m_idx),
        gen_flanking_msa(input_tsv, m_idx),
        gen_flank_lt_msa(input_tsv, m_idx),
        gen_flank_rt_msa(input_tsv, m_idx)
    ]
    return result


def gen_phased_fastas(input_tsv, m1_idx, m2_idx) -> list[str]:
    result: list[str] = [
        # These functions look similar, but are subtly different.
        # I think it is better to keep them separate, to make reasoning within them easier - less cases to keep in mind.
        # Also, their interface is simple - given pandas DataFrame in tsv and module, they return str representing MSA.
        gen_two_good_msa(input_tsv, m1_idx, m2_idx),
        gen_one_good_msa(input_tsv, m1_idx, m2_idx),
        gen_1good_lt_msa(input_tsv, m1_idx, m2_idx),
        gen_1good_rt_msa(input_tsv, m1_idx, m2_idx)
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
    seq_logo = 'true'
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


# --- classes
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


class ChromEnum(enum.Enum):
    X = 'X'
    Y = 'Y'
    NORM = 'NORM'
# --- end classes


def gen_spanning_msa(input_tsv: str, m1: int) -> str:
    df = pd.read_csv(input_tsv, sep='\t')
    n_modules = list(df["n_modules"])[0]

    mask = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Spanning")
    new_df = df[mask]
    if len(new_df) == 0:
        return ""
    raligns = [True] * m1 + [False] * (n_modules - m1)
    realigned_seqs = realign_sequences(new_df, raligns)  # This already contains all information, rest is formatting

    tmp1 = sorted(realigned_seqs, key=lambda x: x[1], reverse=True)
    tmp2 = sorted(tmp1, key=lambda x: x[1][m1].count("-"))
    msa = [(x[0], "-".join(x[1])) for x in tmp2]

    result = msa_to_str(msa)
    return result


def gen_spann_al_msa(input_tsv: str, m1: int, a1: int) -> str:
    df = pd.read_csv(input_tsv, sep='\t')
    msa = _gen_msa_spanning_allele(df, m1, a1)
    result = msa_to_str(msa)
    return result


def gen_flank_lt_msa(input_tsv: str, m1: int) -> str:
    df = pd.read_csv(input_tsv, sep='\t')
    msa = _gen_msa_flanking_left(df, m1)
    result = msa_to_str(msa)
    return result


def gen_flank_rt_msa(input_tsv: str, m1: int) -> str:
    df = pd.read_csv(input_tsv, sep='\t')
    msa = _gen_msa_flanking_right(df, m1)
    result = msa_to_str(msa)
    return result


def gen_flanking_msa(input_tsv: str, m1: int) -> str:
    df = pd.read_csv(input_tsv, sep='\t')
    n_modules = list(df["n_modules"])[0]

    tmp1 = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Flanking")
    tmp2 = (df["module_classes"].apply(lambda x: x.split(",")[m1-1]) == "Missing")
    mask_left = tmp1 & tmp2
    new_df_left = df[mask_left]
    max_left = get_max_widths(new_df_left, n_modules)

    tmp3 = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Flanking")
    tmp4 = (df["module_classes"].apply(lambda x: x.split(",")[m1+1]) == "Missing")
    mask_right = tmp3 & tmp4
    new_df_right = df[mask_right]
    max_right = get_max_widths(new_df_right, n_modules)

    maximums = [max(x, y) for x, y in zip(max_left, max_right)]
    maximums[m1] = max_left[m1] + max_right[m1]

    raligns = [True] * (m1 + 1) + [False] * (n_modules - m1 - 1)
    realigned_seqs_left = []
    for i, row in new_df_left.iterrows():
        aligned_parts = get_aligned_parts(row, raligns, maximums, n_modules)
        item = (row["read_id"], aligned_parts)
        realigned_seqs_left.append(item)
    tmp11 = sorted(realigned_seqs_left, key=lambda x: x[1], reverse=True)
    tmp21 = sorted(tmp11, key=lambda x: x[1][m1].count("-"))
    aln_r = [(x[0], "-".join(x[1])) for x in tmp21]

    raligns = [True] * m1 + [False] * (n_modules - m1)
    realigned_seqs_right = []
    for i, row in new_df_right.iterrows():
        aligned_parts = get_aligned_parts(row, raligns, maximums, n_modules)
        item = (row["read_id"], aligned_parts)
        realigned_seqs_right.append(item)
    tmp31 = sorted(realigned_seqs_right, key=lambda x: x[1], reverse=True)
    tmp41 = sorted(tmp31, key=lambda x: x[1][m1].count("-"))
    aln_l = [(x[0], "-".join(x[1])) for x in tmp41]

    msa = aln_l + [("", "")] + aln_r
    result = msa_to_str(msa)
    return result


def gen_two_good_msa(input_tsv: str, m1: int, m2: int) -> str:
    df = pd.read_csv(input_tsv, sep='\t')
    msa = _gen_msa_two_good(df, m1, m2)
    result = msa_to_str(msa)
    return result


def gen_1good_lt_msa(input_tsv: str, m1: int, m2: int) -> str:
    df = pd.read_csv(input_tsv, sep='\t')
    msa = _gen_msa_one_good_left(df, m1, m2)
    result = msa_to_str(msa)
    return result


def gen_1good_rt_msa(input_tsv: str, m1: int, m2: int) -> str:
    df = pd.read_csv(input_tsv, sep='\t')
    msa = _gen_msa_one_good_right(df, m1, m2)
    result = msa_to_str(msa)
    return result


def gen_one_good_msa(input_tsv: str, m1: int, m2: int) -> str:
    df = pd.read_csv(input_tsv, sep='\t')
    n_modules = list(df["n_modules"])[0]

    tmp1 = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Spanning")
    tmp2 = (df["module_classes"].apply(lambda x: x.split(",")[m2]) == "Flanking")
    mask_left = tmp1 & tmp2
    new_df_left = df[mask_left]
    max_left = get_max_widths(new_df_left, n_modules)

    tmp3 = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Flanking")
    tmp4 = (df["module_classes"].apply(lambda x: x.split(",")[m2]) == "Spanning")
    mask_right = tmp3 & tmp4
    new_df_right = df[mask_right]
    max_right = get_max_widths(new_df_right, n_modules)

    maximums = [max(x, y) for x, y in zip(max_left, max_right)]
    raligns = [True] * (m1 + 1) + [False] * (n_modules - m1 - 1)

    realigned_seqs_right = []
    for i, row in new_df_right.iterrows():
        aligned_parts = get_aligned_parts(row, raligns, maximums, n_modules)
        item = (row["read_id"], aligned_parts)
        realigned_seqs_right.append(item)
    tmp31 = sorted(realigned_seqs_right, key=lambda x: x[1], reverse=True)
    tmp41 = sorted(tmp31, key=lambda x: x[1][m1].count("-"))
    aln_l = [(x[0], "-".join(x[1])) for x in tmp41]

    realigned_seqs_left = []
    for i, row in new_df_left.iterrows():
        aligned_parts = get_aligned_parts(row, raligns, maximums, n_modules)
        item = (row["read_id"], aligned_parts)
        realigned_seqs_left.append(item)
    tmp11 = sorted(realigned_seqs_left, key=lambda x: x[1], reverse=True)
    tmp21 = sorted(tmp11, key=lambda x: x[1][m2].count("-"))
    aln_r = [(x[0], "-".join(x[1])) for x in tmp21]

    msa = aln_l + [("", "")] + aln_r
    result = msa_to_str(msa)
    return result


def msa_to_str(msa: MSA) -> str:
    string_tmp = []
    for annot_name, align in msa:
        string_tmp.append(f">{annot_name}\n{align}\n")
    string = "".join(string_tmp)
    return string


def _gen_msa_spanning_allele(df: pd.DataFrame, m1: int, a1: int) -> MSA:
    n_modules = list(df["n_modules"])[0]
    tmp1 = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Spanning")
    tmp2 = (df["module_repetitions"].apply(lambda x: int(x.split(",")[m1])) == a1)
    mask = tmp1 & tmp2
    new_df = df[mask]
    if len(new_df) == 0:
        return []
    raligns = [True] * m1 + [False] * (n_modules - m1)
    realigned_seqs = realign_sequences(new_df, raligns)  # This already contains all information, rest is formatting

    tmp11 = sorted(realigned_seqs, key=lambda x: x[1], reverse=True)
    tmp21 = sorted(tmp11, key=lambda x: x[1][m1].count("-"))
    result = [(x[0], "-".join(x[1])) for x in tmp21]
    return result


def _gen_msa_flanking_left(df: pd.DataFrame, m1: int) -> MSA:
    n_modules = list(df["n_modules"])[0]
    tmp1 = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Flanking")
    tmp2 = (df["module_classes"].apply(lambda x: x.split(",")[m1+1]) == "Missing")
    mask = tmp1 & tmp2
    new_df = df[mask]
    if len(new_df) == 0:
        return []
    raligns = [True] * m1 + [False] * (n_modules - m1)
    realigned_seqs = realign_sequences(new_df, raligns)

    tmp1 = sorted(realigned_seqs, key=lambda x: x[1], reverse=True)
    tmp2 = sorted(tmp1, key=lambda x: x[1][m1].count("-"))
    result = [(x[0], "-".join(x[1])) for x in tmp2]
    return result


def _gen_msa_flanking_right(df: pd.DataFrame, m1: int) -> MSA:
    n_modules = list(df["n_modules"])[0]
    tmp1 = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Flanking")
    tmp2 = (df["module_classes"].apply(lambda x: x.split(",")[m1-1]) == "Missing")
    mask = tmp1 & tmp2
    new_df = df[mask]
    if len(new_df) == 0:
        return []
    raligns = [True] * (m1 + 1) + [False] * (n_modules - m1 - 1)
    realigned_seqs = realign_sequences(new_df, raligns)

    tmp1 = sorted(realigned_seqs, key=lambda x: x[1], reverse=True)
    tmp2 = sorted(tmp1, key=lambda x: x[1][m1].count("-"))
    result = [(x[0], "-".join(x[1])) for x in tmp2]
    return result


def _gen_msa_flanking(df: pd.DataFrame, m1: int) -> MSA:
    aln_r = _gen_msa_flanking_right(df, m1)
    aln_l = _gen_msa_flanking_left(df, m1)
    return aln_r + [("", "")] + aln_l


def _gen_msa_two_good(df: pd.DataFrame, m1: int, m2: int) -> MSA:
    n_modules = list(df["n_modules"])[0]
    tmp1 = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Spanning")
    tmp2 = (df["module_classes"].apply(lambda x: x.split(",")[m2]) == "Spanning")
    mask = tmp1 & tmp2
    new_df = df[mask]
    if len(new_df) == 0:
        return []
    raligns = [True] * m1 + [False] * (n_modules - m1)
    realigned_seqs = realign_sequences(new_df, raligns)

    tmp1 = sorted(realigned_seqs, key=lambda x: x[1], reverse=True)
    tmp2 = sorted(tmp1, key=lambda x: x[1][m2].count("-"))
    tmp3 = sorted(tmp2, key=lambda x: x[1][m1].count("-"))
    result = [(x[0], "-".join(x[1])) for x in tmp3]
    return result


def _gen_msa_one_good_left(df: pd.DataFrame, m1: int, m2: int) -> MSA:
    n_modules = list(df["n_modules"])[0]
    tmp1 = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Spanning")
    tmp2 = (df["module_classes"].apply(lambda x: x.split(",")[m2]) == "Flanking")
    mask = tmp1 & tmp2
    new_df = df[mask]
    if len(new_df) == 0:
        return []
    raligns = [True] * m1 + [False] * (n_modules - m1)
    realigned_seqs = realign_sequences(new_df, raligns)  # This already contains all information, rest is formatting

    tmp1 = sorted(realigned_seqs, key=lambda x: x[1], reverse=True)
    tmp2 = sorted(tmp1, key=lambda x: x[1][m1].count("-"))
    result = [(x[0], "-".join(x[1])) for x in tmp2]
    return result


def _gen_msa_one_good_right(df: pd.DataFrame, m1: int, m2: int) -> MSA:
    n_modules = list(df["n_modules"])[0]
    tmp1 = (df["module_classes"].apply(lambda x: x.split(",")[m1]) == "Flanking")
    tmp2 = (df["module_classes"].apply(lambda x: x.split(",")[m2]) == "Spanning")
    mask = tmp1 & tmp2
    new_df = df[mask]
    if len(new_df) == 0:
        return []
    raligns = [True] * m2 + [False] * (n_modules - m2)
    realigned_seqs = realign_sequences(new_df, raligns)  # This already contains all information, rest is formatting

    tmp1 = sorted(realigned_seqs, key=lambda x: x[1], reverse=True)
    tmp2 = sorted(tmp1, key=lambda x: x[1][m2].count("-"))
    result = [(x[0], "-".join(x[1])) for x in tmp2]
    return result


def realign_sequences(df: pd.DataFrame, raligns: list[bool]):
    n = len(raligns)
    maximums = get_max_widths(df, n)
    result = []
    for i, row in df.iterrows():
        aligned_parts = get_aligned_parts(row, raligns, maximums, n)
        item = (row["read_id"], aligned_parts)
        result.append(item)

    return result


def get_aligned_parts(row: pd.Series, raligns: list[bool], max_widths: list[int], n: int) -> list[str]:
    aligned_parts = []
    parts = row["module_sequences"].split(",")
    for j in range(n):
        if raligns[j]:
            x = "{0:->{w}}".format(parts[j], w=max_widths[j])
        else:
            x = "{0:-<{w}}".format(parts[j], w=max_widths[j])
        aligned_parts.append(x)
    return aligned_parts


def get_max_widths(df: pd.DataFrame, n: int) -> list[int]:
    maximums = [0] * n
    for i, row in df.iterrows():
        parts = row["module_sequences"].split(",")
        for j in range(n):
            maximums[j] = max(len(parts[j]), maximums[j])
    return maximums


def highlight_subpart(seq: str, highlight: int | list[int]) -> tuple[str, str]:
    if highlight is None:
        return seq, ''

    str_part = []
    highlight1 = np.array(highlight)
    split = [f'{s}]' for s in seq.split(']') if s != '']
    for h in highlight1:
        str_part.append(split[h])
        split[h] = f'<b><u>{split[h]}</u></b>'
    return ''.join(split), ''.join(str_part)


def create_motif(df: pd.DataFrame, is_male: bool) -> Motif:
    motif_str = df[MOTIF_COLUMN_NAME].iloc[0]
    name = None if 'name' not in df.columns or df.iloc[0]['name'] in ['None', ''] else df.iloc[0]['name']
    motif_class = Motif(motif_str, name)
    motif_class.monoallelic = is_male and chrom_from_string(motif_class.chrom) in [ChromEnum.X, ChromEnum.Y]
    return motif_class


def chrom_from_string(chrom_str: str) -> ChromEnum:
    return (
        ChromEnum.X if chrom_str in ['chrX', 'NC_000023'] else
        ChromEnum.Y if chrom_str in ['chrY', 'NC_000024'] else
        ChromEnum.NORM
    )


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
