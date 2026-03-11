from __future__ import annotations

import re
from enum import Enum

from src_new.constants import MOTIF_COLUMN_NAME, MOTIF_COLUMN_ID


def generate_motif_stats(motif_table, male):
    motif = Motif(motif_table[MOTIF_COLUMN_NAME].iloc[0], motif_table[MOTIF_COLUMN_ID].iloc[0], male)
    return motif.get_motif_stats()


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

        motif_monoallelic = male and ChromEnum.from_string(chrom) in [ChromEnum.X, ChromEnum.Y]

        # store members
        self.nomenclature: str = nomenclature
        self.name: str = name
        self.chrom: str = chrom
        self.start: int = int(start)
        self.end: int = int(end)
        self.modules: list[tuple[str, int]] = modules
        self.monoallelic: bool = motif_monoallelic

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

    def get_motif_stats(self) -> dict:
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "modules": self.modules
        }

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


class ChromEnum(Enum):
    X = 'X'
    Y = 'Y'
    NORM = 'NORM'

    @staticmethod
    def from_string(chrom_str: str) -> ChromEnum:
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
