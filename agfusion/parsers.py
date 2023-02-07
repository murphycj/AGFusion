"""
Parses output files from fusion-finding algorithms
"""

import csv
import re
import sys

import pandas as pd


class _Parser:
    """Base parser."""

    def __init__(self, logger):
        self.fusions = []
        self.iterator = 0
        self.logger = logger

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterator >= len(self.fusions):
            raise StopIteration

        self.iterator += 1
        return self.fusions[self.iterator - 1]

    def _check_data(self):
        if len(self.fusions) == 0:
            self.logger.error("Read 0 fusions from the file! Exiting...")
        else:
            self.logger.info(f"Read {len(self.fusions)} fusions from the file.")

    next = __next__


class STARFusion(_Parser):
    """STARFusion parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        with open(infile, "r") as fin:
            reader = csv.DictReader(fin, delimiter="\t")
            if not ("#FusionName" in reader.fieldnames or "#fusion_name" in reader.fieldnames):
                raise AssertionError(
                    "Unrecognized STAR-Fusion input for first column "
                    + "in header. Should be #FusionName or #fusion_name."
                )

            assert "LeftGene" in reader.fieldnames, "Unrecognized STAR-Fusion input"
            assert "LeftBreakpoint" in reader.fieldnames, "Unrecognized " + "STAR-Fusion input"
            assert "RightGene" in reader.fieldnames, "Unrecognized STAR-Fusion input"
            assert "RightBreakpoint" in reader.fieldnames, "Unrecognized " + "STAR-Fusion input"

            for line in reader:
                gene_5prime = line["LeftGene"].split("^")[1].split(".")[0]
                gene_5prime_name = line["LeftGene"].split("^")[0]
                gene_5prime_junction = int(line["LeftBreakpoint"].split(":")[1])
                gene_3prime = line["RightGene"].split("^")[1].split(".")[0]
                gene_3prime_name = line["RightGene"].split("^")[0]
                gene_3prime_junction = int(line["RightBreakpoint"].split(":")[1])
                self.fusions.append(
                    {
                        "gene5prime": gene_5prime,
                        "gene3prime": gene_3prime,
                        "alternative_name_5prime": gene_5prime_name,
                        "alternative_name_3prime": gene_3prime_name,
                        "gene5prime_junction": gene_5prime_junction,
                        "gene3prime_junction": gene_3prime_junction,
                    }
                )

        self._check_data()


class Arriba(_Parser):
    """Arriba parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        data = pd.read_csv(infile, delimiter="\t")
        data.columns = [i.replace("#", "") for i in data.columns]

        cols = ["gene1", "gene_id1", "gene2", "gene_id2", "breakpoint1", "breakpoint2"]

        assert all(i in data.columns for i in cols), "Unrecognized Arriba input"

        for i in data.index:
            if not data.at[i, "gene1"]:
                continue

            gene_5prime = data.at[i, "gene_id1"].split(".")[0]
            gene_5prime_name = data.at[i, "gene1"]
            gene_5prime_junction = int(data.at[i, "breakpoint1"].split(":")[1])
            gene_3prime = data.at[i, "gene_id2"].split(".")[0]
            gene_3prime_name = data.at[i, "gene2"]
            gene_3prime_junction = int(data.at[i, "breakpoint2"].split(":")[1])
            self.fusions.append(
                {
                    "gene5prime": gene_5prime,
                    "gene3prime": gene_3prime,
                    "alternative_name_5prime": gene_5prime_name,
                    "alternative_name_3prime": gene_3prime_name,
                    "gene5prime_junction": gene_5prime_junction,
                    "gene3prime_junction": gene_3prime_junction,
                }
            )

        self._check_data()


class EricScript(_Parser):
    """EricScript parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        fin = open(infile, "r")
        for line in fin.readlines():
            if re.findall(r"^GeneName1", line):
                line = line.strip().split("\t")
                for i, j in zip(
                    [3, 6, 8, 9],
                    ["Breakpoint1", "Breakpoint2", "EnsemblGene1", "EnsemblGene2"],
                ):
                    assert line[i] == j, f"Unrecognized EricScript input: {j}"
            else:
                line = line.strip().split("\t")

                gene_5prime_name = line[8]
                gene_5prime_junction = int(line[3])
                gene_3prime_name = line[9]
                gene_3prime_junction = int(line[6])
                self.fusions.append(
                    {
                        "gene5prime": None,
                        "gene3prime": None,
                        "alternative_name_5prime": gene_5prime_name,
                        "alternative_name_3prime": gene_3prime_name,
                        "gene5prime_junction": gene_5prime_junction,
                        "gene3prime_junction": gene_3prime_junction,
                    }
                )
        fin.close()

        self._check_data()


class FusionCatcher(_Parser):
    """FusionCatcher parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        fin = open(infile, "r")
        for line in fin.readlines():
            if re.findall(r"^Gene_1_symbol", line):
                line = line.rstrip().split("\t")
                assert (
                    line[8] == "Fusion_point_for_gene_1(5end_fusion_partner)"
                ), "Unrecognized FusionCatcher input"
                assert (
                    line[9] == "Fusion_point_for_gene_2(3end_fusion_partner)"
                ), "Unrecognized FusionCatcher input"
                assert (
                    line[10] == "Gene_1_id(5end_fusion_partner)"
                ), "Unrecognized FusionCatcher input"
                assert (
                    line[11] == "Gene_2_id(3end_fusion_partner)"
                ), "Unrecognized FusionCatcher input"
                continue

            line = line.strip().split("\t")
            self.fusions.append(
                {
                    "gene5prime": line[10],
                    "gene3prime": line[11],
                    "alternative_name_5prime": line[0],
                    "alternative_name_3prime": line[1],
                    "gene5prime_junction": int(line[8].split(":")[1]),
                    "gene3prime_junction": int(line[9].split(":")[1]),
                }
            )
        fin.close()

        self._check_data()


class FusionHunter(_Parser):
    """FusionHunter parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        gene1 = gene2 = None
        gene1_junction = gene2_junction = None

        fin = open(infile, "r")
        for line in fin.readlines():

            if re.findall(r"^# Fusion:", line):

                if gene1 is not None and gene2 is not None:
                    self.fusions.append(
                        {
                            "gene5prime": None,
                            "gene3prime": None,
                            "alternative_name_5prime": gene1,
                            "alternative_name_3prime": gene2,
                            "gene5prime_junction": int(gene1_junction),
                            "gene3prime_junction": int(gene2_junction),
                        }
                    )

                strands = re.findall(r"(?<=\[)(.*)(?=\])", line)
                assert (
                    len(strands) == 1
                ), "Unrecognized FusionHunter input. Incorrect strand information."
                gene1_strand = strands[0][0]
                gene2_strand = strands[0][1]

            elif re.findall(r"^--", line):
                # new breakpoint
                if gene1 is not None and gene2 is not None:
                    self.fusions.append(
                        {
                            "gene5prime": None,
                            "gene3prime": None,
                            "alternative_name_5prime": gene1,
                            "alternative_name_3prime": gene2,
                            "gene5prime_junction": int(gene1_junction),
                            "gene3prime_junction": int(gene2_junction),
                        }
                    )

            elif re.findall(r"^->", line):
                junctions = re.findall(r"(chr[0-9]*):([0-9]*)-([0-9]*)", line)
                assert len(junctions) == 2, (
                    "Unrecognized FusionHunter " + "input. Incorrect junction information."
                )

                gene1 = gene2 = None
                gene1_junction = gene2_junction = None

                gene1, gene2 = re.findall(r"[A-Z0-9]*\sx\s[A-Z0-9a-z]*(?=\t)", line)[0].split(" x ")
                assert (
                    gene1 is not None and gene2 is not None
                ), "Unrecognized FusionHunter input. Incorrect gnee information."

                if gene1_strand == "+":
                    gene1_junction = junctions[0][2]
                else:
                    gene1_junction = junctions[0][1]

                if gene2_strand == "+":
                    gene2_junction = junctions[1][1]
                else:
                    gene2_junction = junctions[1][2]
        fin.close()

        if gene1 is not None and gene2 is not None:
            self.fusions.append(
                {
                    "gene5prime": None,
                    "gene3prime": None,
                    "alternative_name_5prime": gene1,
                    "alternative_name_3prime": gene2,
                    "gene5prime_junction": int(gene1_junction),
                    "gene3prime_junction": int(gene2_junction),
                }
            )

        self._check_data()


class FusionMap(_Parser):
    """FusionMap parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        fin = open(infile, "r")
        for line in fin.readlines():
            if re.findall(r"^FusionID", line):
                line = line.strip().split("\t")
                for i, j in zip(
                    [6, 8, 9, 13],
                    ["Position1", "Position2", "KnownGene1", "KnownGene2"],
                ):
                    assert (
                        line[i] == j
                    ), f"Unrecognized FusionMap input {j} not in expected position."
            else:
                line = line.strip().split("\t")

                gene_5prime_name = line[9]
                gene_5prime_junction = int(line[6])
                gene_3prime_name = line[13]
                gene_3prime_junction = int(line[8])
                self.fusions.append(
                    {
                        "gene5prime": None,
                        "gene3prime": None,
                        "alternative_name_5prime": gene_5prime_name,
                        "alternative_name_3prime": gene_3prime_name,
                        "gene5prime_junction": gene_5prime_junction,
                        "gene3prime_junction": gene_3prime_junction,
                    }
                )
        fin.close()

        self._check_data()


class MapSplice(_Parser):
    """MapSplice parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        fin = open(infile, "r")
        for line in fin.readlines():
            if re.findall(r"^chrom", line):
                line = line.strip().split("\t")
                for i, j in zip(
                    [1, 2, 60, 61],
                    [
                        "doner_end",
                        "acceptor_start",
                        "annotated_gene_donor",
                        "annotated_gene_acceptor",
                    ],
                ):
                    assert (
                        line[i] == j
                    ), f"Unrecognized MapSplice input. {j} header not in expected position."
            else:
                line = line.strip().split("\t")

                gene_5prime_name = line[60]
                gene_5prime_junction = int(line[1])
                gene_3prime_name = line[61]
                gene_3prime_junction = int(line[2])
                self.fusions.append(
                    {
                        "gene5prime": None,
                        "gene3prime": None,
                        "alternative_name_5prime": gene_5prime_name,
                        "alternative_name_3prime": gene_3prime_name,
                        "gene5prime_junction": gene_5prime_junction,
                        "gene3prime_junction": gene_3prime_junction,
                    }
                )
        fin.close()

        self._check_data()


class TopHatFusion(_Parser):
    """TopHatFusion parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        fin = open(infile, "r")
        for line in fin.readlines():
            line = line.strip().split("\t")

            gene_5prime = line[1]
            gene_5prime_name = line[1]
            gene_5prime_junction = int(line[3])
            gene_3prime = line[4]
            gene_3prime_name = line[4]
            gene_3prime_junction = int(line[6])
            self.fusions.append(
                {
                    "gene5prime": gene_5prime,
                    "gene3prime": gene_3prime,
                    "alternative_name_5prime": gene_5prime_name,
                    "alternative_name_3prime": gene_3prime_name,
                    "gene5prime_junction": gene_5prime_junction,
                    "gene3prime_junction": gene_3prime_junction,
                }
            )
        fin.close()

        self._check_data()


class DeFuse(_Parser):
    """DeFuse parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        data_indices = {
            "gene5prime": None,
            "gene3prime": None,
            "gene5prime_junction": None,
            "gene3prime_junction": None,
        }

        fin = open(infile, "r")
        for line in fin.readlines():
            if re.findall(r"^cluster_id", line):
                line = line.strip().split("\t")
                for column in data_indices:
                    try:
                        data_indices[column] = line.index(column)
                    except ValueError:
                        logger.error(
                            f"Unrecognized {self.__class__.__name__} input! "
                            f"Cannot find {column} column."
                        )
                        sys.exit(1)
                continue
            if not any(i is None for i in data_indices):
                line = line.strip().split("\t")
                self.fusions.append(
                    {
                        "gene5prime": line[data_indices["gene5prime"]],
                        "gene3prime": line[data_indices["gene3prime"]],
                        "gene5prime_junction": int(line[data_indices["gene5prime_junction"]]),
                        "gene3prime_junction": int(line[data_indices["gene3prime_junction"]]),
                    }
                )
        fin.close()

        self._check_data()


class Chimerascan(_Parser):
    """Chimerascan parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        data_indices = {
            "genes5p": 12,
            "genes3p": 13,
            "start5p": 1,
            "end5p": 2,
            "start3p": 4,
            "end3p": 5,
            "strand5p": 8,
            "strand3p": 9,
        }

        fin = open(infile, "r")
        for line in fin.readlines():
            if re.findall(r"#chrom5p", line):
                line = line.strip().split("\t")
                for column in data_indices:
                    try:
                        data_indices[column] = line.index(column)
                    except ValueError:
                        logger.error(
                            f"Unrecognized {self.__class__.__name__} input! "
                            f"Cannot find {column} column."
                        )
                        sys.exit(1)
                continue

            line = line.strip().split("\t")

            if line[data_indices["strand5p"]] == "+":
                gene1_junction = line[data_indices["end5p"]]
            else:
                gene1_junction = line[data_indices["start5p"]]

            if line[data_indices["strand3p"]] == "+":
                gene2_junction = line[data_indices["start3p"]]
            else:
                gene2_junction = line[data_indices["end3p"]]

            self.fusions.append(
                {
                    "gene5prime": line[data_indices["genes5p"]].split(","),
                    "gene3prime": line[data_indices["genes3p"]].split(","),
                    "gene5prime_junction": int(gene1_junction),
                    "gene3prime_junction": int(gene2_junction),
                }
            )
        fin.close()

        self._check_data()


class ChimeRScope(_Parser):
    """ChimeRScope parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        data_indices = {
            "Gene1": 2,
            "Gene2": 4,
            "Gene1_fusionPoint": 7,
            "Gene2_fusionPoint": 9,
        }

        fin = open(infile, "r")
        for line in fin.readlines():
            if re.findall(r"^ConfidentScore", line):
                line = line.strip().split("\t")
                for column in data_indices:
                    try:
                        data_indices[column] = line.index(column)
                    except ValueError:
                        logger.error(
                            f"Unrecognized {self.__class__.__name__} input! "
                            f"Cannot find {column} column."
                        )
                        sys.exit(1)
                continue

            line = line.strip().split("\t")
            self.fusions.append(
                {
                    "gene5prime": line[data_indices["Gene1"]],
                    "gene3prime": line[data_indices["Gene2"]],
                    "gene5prime_junction": line[data_indices["Gene1_fusionPoint"]],
                    "gene3prime_junction": line[data_indices["Gene2_fusionPoint"]],
                }
            )
        fin.close()

        self._check_data()


class JAFFA(_Parser):
    """JAFFA parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        data_indices = {"base1": 7, "base2": 9, "fusion genes": 1}

        fin = open(infile, "r")
        for line in fin.readlines():
            if re.findall(r"sample", line):
                line = line.strip().replace('"', "").split(",")
                for column in data_indices:
                    try:
                        data_indices[column] = line.index(column)
                    except ValueError:
                        logger.error(
                            f"Unrecognized {self.__class__.__name__} input! "
                            f"Cannot find {column} column."
                        )
                        sys.exit(1)
                continue

            line = line.strip().replace('"', "").split(",")
            self.fusions.append(
                {
                    "gene5prime": line[data_indices["fusion genes"]].split(":")[0],
                    "gene3prime": line[data_indices["fusion genes"]].split(":")[1],
                    "gene5prime_junction": int(line[data_indices["base1"]]),
                    "gene3prime_junction": int(line[data_indices["base2"]]),
                }
            )
        fin.close()

        self._check_data()


class Bellerophontes(_Parser):
    """Bellerophontes parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        fin = open(infile, "r")
        for line in fin.readlines():

            line = line.strip().split("\t")
            if len(line) <= 2:
                continue

            self.fusions.append(
                {
                    "gene5prime": line[0],
                    "gene3prime": line[4],
                    "gene5prime_junction": int(line[9]),
                    "gene3prime_junction": int(line[11]),
                }
            )
        fin.close()

        self._check_data()


class BreakFusion(_Parser):
    """BreakFusion parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        data_indices = {"POS1": 1, "POS2": 4, "RefseqGene": 11}

        fin = open(infile, "r")
        for line in fin.readlines():
            if re.findall(r"CHR1", line):
                line = line.strip().split("\t")
                for column in data_indices:
                    try:
                        data_indices[column] = line.index(column)
                    except ValueError:
                        logger.error(
                            f"Unrecognized {self.__class__.__name__} input! "
                            f"Cannot find {column} column."
                        )
                        sys.exit(1)
                continue

            if re.findall(r"Fusion", line):
                line = line.strip().split("\t")
                genes = line[data_indices["RefseqGene"]]
                genes = re.findall(r"(?<=Gene:)(.*?)(?=,)", genes)
                if len(genes) != 1:
                    line_out = "\t".join(line)
                    self.logger.error(
                        f"Could not parse genes from BreakFusion input line: {line_out}"
                    )
                genes = genes[0].split("|")
                if len(genes) != 2:
                    line_out = "\t".join(line)
                    self.logger.error(
                        f"Could not find two genes for BreakFusion input line: {line_out}"
                    )

                self.fusions.append(
                    {
                        "gene5prime": genes[0],
                        "gene3prime": genes[1],
                        "gene5prime_junction": int(line[data_indices["POS1"]]),
                        "gene3prime_junction": int(line[data_indices["POS2"]]),
                    }
                )
        fin.close()

        self._check_data()


class InFusion(_Parser):
    """InFusion parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        data_indices = {"break_pos1": 2, "break_pos2": 5, "genes_1": 9, "genes_2": 10}

        fin = open(infile, "r")
        n = 0
        for line in fin.readlines():
            n += 1
            if re.findall(r"#id", line):
                line = line.strip().split("\t")
                for column in data_indices:
                    try:
                        data_indices[column] = line.index(column)
                    except ValueError:
                        logger.error(
                            f"Unrecognized {self.__class__.__name__} input! "
                            f"Cannot find {column} column."
                        )
                        sys.exit(1)
                continue

            line = line.strip().split("\t")

            if line[data_indices["genes_1"]] == "none" or line[data_indices["genes_2"]] == "none":
                self.logger.warn(
                    f"Skipping fusion on line {n} because one or more "
                    + "of the provided gene names under 'gene_1' and"
                    + " 'gene_2' is listed as 'none'."
                )
                continue

            self.fusions.append(
                {
                    "gene5prime": line[data_indices["genes_1"]].split(";"),
                    "gene3prime": line[data_indices["genes_2"]].split(";"),
                    "gene5prime_junction": int(line[data_indices["break_pos1"]]),
                    "gene3prime_junction": int(line[data_indices["break_pos2"]]),
                }
            )
        fin.close()

        self._check_data()


class FusionInspector(_Parser):
    """FusionInspector parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        data = pd.read_csv(infile, delimiter="\t")
        data.columns = [i.replace("#", "") for i in data.columns]

        cols = ["LeftGene", "LeftBreakpoint", "RightGene", "RightBreakpoint"]
        assert all(
            i in data.columns for i in cols
        ), "Unrecognized FusionInspector input. Could not find all columns: " + ",".join(cols)

        for i in data.index:

            gene_5prime = data.at[i, "LeftGene"].split("^")[1].split(".")[0]
            gene_5prime_name = data.at[i, "LeftGene"].split("^")[0]
            gene_5prime_junction = int(data.at[i, "LeftBreakpoint"].split(":")[1])
            gene_3prime = data.at[i, "RightGene"].split("^")[1].split(".")[0]
            gene_3prime_name = data.at[i, "RightGene"].split("^")[0]
            gene_3prime_junction = int(data.at[i, "RightBreakpoint"].split(":")[1])
            self.fusions.append(
                {
                    "gene5prime": gene_5prime,
                    "gene3prime": gene_3prime,
                    "alternative_name_5prime": gene_5prime_name,
                    "alternative_name_3prime": gene_3prime_name,
                    "gene5prime_junction": gene_5prime_junction,
                    "gene3prime_junction": gene_3prime_junction,
                }
            )

        self._check_data()


class LongGf(_Parser):
    """LongGf parser."""

    def __init__(self, infile, logger):
        super().__init__(logger)

        fin = open(infile, "r")
        for line in fin.readlines():
            if re.findall(r"^GF", line):
                line = line.rstrip().replace("\t", " ").split(" ")

                gene_5prime, gene_3prime = line[1].split(":")

                assert gene_5prime and gene_3prime, "Could not parse LongGf genes: " + line[1]

                gene_5prime_chrom, gene_5prime_junction = line[8].split(":")

                assert gene_5prime_chrom and gene_5prime_junction, (
                    "Could not parse LongGf junction: " + line[8]
                )

                gene_3prime_chrom, gene_3prime_junction = line[10].split(":")

                assert gene_3prime_chrom and gene_3prime_junction, (
                    "Could not parse LongGf junction: " + line[10]
                )

                self.fusions.append(
                    {
                        "gene5prime": gene_5prime,
                        "gene3prime": gene_3prime,
                        "alternative_name_5prime": gene_5prime,
                        "alternative_name_3prime": gene_3prime,
                        "gene5prime_junction": int(gene_5prime_junction),
                        "gene3prime_junction": int(gene_3prime_junction),
                    }
                )
        fin.close()

        self._check_data()


parsers = {
    "arriba": Arriba,
    "bellerophontes": Bellerophontes,
    "breakfusion": BreakFusion,
    "chimerascan": Chimerascan,
    "chimerscope": ChimeRScope,
    "defuse": DeFuse,
    "ericscript": EricScript,
    "fusioncatcher": FusionCatcher,
    "fusionhunter": FusionHunter,
    "fusionmap": FusionMap,
    "fusioninspector": FusionInspector,
    "infusion": InFusion,
    "jaffa": JAFFA,
    "longgf": LongGf,
    "starfusion": STARFusion,
    "tophatfusion": TopHatFusion,
    "mapsplice": MapSplice,
    # 'fusionseq':FusionSeq,
    # 'prada':Prada,
    # 'gfusion':GFusion,
    # 'machete':Machete,
    # 'nfuse':nFuse,
    # 'soapfuse':SOAPfuse,
}
