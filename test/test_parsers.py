"""Unit tests.
"""

import unittest
from os.path import abspath, curdir, join

import pyensembl

from agfusion import database, model, parsers

data_mouse = pyensembl.EnsemblRelease(84, "mouse")
db_mouse = database.AGFusionDB(abspath(join(curdir, "agfusion.mus_musculus.84.db")))
db_mouse.build = "mus_musculus_84"

data_human = pyensembl.EnsemblRelease(75, "human")
db_human = database.AGFusionDB(abspath(join(curdir, "agfusion.homo_sapiens.75.db")))
db_human.build = "homo_sapiens_75"

data_human_hg38 = pyensembl.EnsemblRelease(111, "human")
db_human_hg38 = database.AGFusionDB(abspath(join(curdir, "agfusion.homo_sapiens.111.db")))
db_human_hg38.build = "homo_sapiens_111"


BASEDIR = "./data/FusionsFindingAlgorithms"


class TestFusionCatcher(unittest.TestCase):
    """Test parse FusionCatcher parse."""

    def test_parse(self):
        """Test basic parsing."""

        all_fusions = [
            "Adamts9_Ano2",
            "Trp53_Sat2",
            "1700112E06Rik_Runx1",
            "Runx1_1700112E06Rik",
            "Rell1_Lhfpl3",
            "Phc1_Smarca2",
            "Lrrc8d_Gbp11",
            "C920009B18Rik_H60b",
        ]
        for fusion in parsers.parsers["fusioncatcher"](
            f"{BASEDIR}/FusionCatcher/final-list_candidate-fusion-genes.txt",
            db_mouse.logger,
        ):
            fusion = model.Fusion(
                gene5prime=fusion["gene5prime"],
                gene5primejunction=fusion["gene5prime_junction"],
                gene3prime=fusion["gene3prime"],
                gene3primejunction=fusion["gene3prime_junction"],
                db=db_mouse,
                pyensembl_data=data_mouse,
                protein_databases=["pfam"],
                noncanonical=False,
            )
            assert fusion.name in all_fusions, f"{fusion.name} not in list!"


class TestArriba(unittest.TestCase):
    """Test parse Arriba parse."""

    def test_parse(self):
        """Test basic parsing."""

        all_fusions = [
            "BCR_ABL1",
        ]
        for fusion in parsers.parsers["arriba"](
            f"{BASEDIR}/Arriba/fusions.tsv",
            db_human.logger,
        ):
            fusion = model.Fusion(
                gene5prime=fusion["gene5prime"],
                gene5primejunction=fusion["gene5prime_junction"],
                gene3prime=fusion["gene3prime"],
                gene3primejunction=fusion["gene3prime_junction"],
                db=db_human,
                pyensembl_data=data_human,
                protein_databases=["pfam"],
                noncanonical=False,
            )
            assert fusion.name in all_fusions, f"{fusion.name} not in list!"


class TestSTARFusion(unittest.TestCase):
    """Test parse STAR Fusion"""

    def test_basic(self):
        """Test basic parsing."""

        all_fusions = ["ACACA_STAC2", "RPS6KB1_SNF8"]
        for fusion in parsers.parsers["starfusion"](
            f"{BASEDIR}/STARFusion/" + "star-fusion.fusion_candidates.final.abridged",
            db_human.logger,
        ):
            fusion = model.Fusion(
                gene5prime=fusion["gene5prime"],
                gene5primejunction=fusion["gene5prime_junction"],
                gene3prime=fusion["gene3prime"],
                gene3primejunction=fusion["gene3prime_junction"],
                db=db_human,
                pyensembl_data=data_human,
                protein_databases=["pfam"],
                noncanonical=False,
            )
            assert fusion.name in all_fusions, f"{fusion.name} not in list!"

    def test_with_coding_effect(self):
        """Test parse output with coding effect."""

        all_fusions = ["ARID3B_MYCNUT", "ARID3B_MYCN", "TVP23C_CDRT4"]
        for fusion in parsers.parsers["starfusion"](
            f"{BASEDIR}/STARFusion/" + "star-fusion.fusion_predictions.abridged.coding_effect.tsv",
            db_human_hg38.logger,
        ):
            fusion = model.Fusion(
                gene5prime=fusion["gene5prime"],
                gene5primejunction=fusion["gene5prime_junction"],
                gene3prime=fusion["gene3prime"],
                gene3primejunction=fusion["gene3prime_junction"],
                db=db_human_hg38,
                pyensembl_data=data_human_hg38,
                protein_databases=["pfam"],
                noncanonical=False,
            )
            assert fusion.name in all_fusions, f"{fusion.name} not in list!"


class TestLongGf(unittest.TestCase):
    """Test parse LongGF"""

    def test_parse_mouse(self):
        """Test basic parsing."""

        all_fusions = ["Mocos_Rprd1a", "Ubc_Ubb", "Ubc_Gm11808", "Gm21887_Gm47283"]
        for fusion in parsers.parsers["longgf"](
            f"{BASEDIR}/LongGF/fusions_mouse.log",
            db_mouse.logger,
        ):
            fusion = model.Fusion(
                gene5prime=fusion["gene5prime"],
                gene5primejunction=fusion["gene5prime_junction"],
                gene3prime=fusion["gene3prime"],
                gene3primejunction=fusion["gene3prime_junction"],
                db=db_mouse,
                pyensembl_data=data_mouse,
                protein_databases=["pfam"],
                noncanonical=False,
            )
            assert fusion.name in all_fusions, f"{fusion.name} not in list!"

    def test_parse_human(self):
        """Test basic parsing."""

        all_fusions = ["BCAS4_BCAS3", "HNRNPC_ACIN1"]
        for fusion in parsers.parsers["longgf"](
            f"{BASEDIR}/LongGF/fusions_hg38.log",
            db_human_hg38.logger,
        ):
            fusion = model.Fusion(
                gene5prime=fusion["gene5prime"],
                gene5primejunction=fusion["gene5prime_junction"],
                gene3prime=fusion["gene3prime"],
                gene3primejunction=fusion["gene3prime_junction"],
                db=db_human_hg38,
                pyensembl_data=data_human_hg38,
                protein_databases=["pfam"],
                noncanonical=False,
            )
            assert fusion.name in all_fusions, f"{fusion.name} not in list!"


class TestFusionInspector(unittest.TestCase):
    """Test parse FusionInspector"""

    def test_parse_human(self):
        """Test basic parsing."""

        all_fusions = ["ENSG00000282885_TPM3", "STAT3_ENSG00000282885"]

        for fusion in parsers.parsers["fusioninspector"](
            f"{BASEDIR}/FusionInspector/test.FusionInspector.fusions.abridged.txt",
            db_human_hg38.logger,
        ):
            fusion = model.Fusion(
                gene5prime=fusion["gene5prime"],
                gene5primejunction=fusion["gene5prime_junction"],
                gene3prime=fusion["gene3prime"],
                gene3primejunction=fusion["gene3prime_junction"],
                db=db_human_hg38,
                pyensembl_data=data_human_hg38,
                protein_databases=["pfam"],
                noncanonical=False,
            )
            assert fusion.name in all_fusions, f"{fusion.name} not in list!"

        for fusion in parsers.parsers["fusioninspector"](
            f"{BASEDIR}/FusionInspector/test.FusionInspector.fusions.txt",
            db_human_hg38.logger,
        ):
            fusion = model.Fusion(
                gene5prime=fusion["gene5prime"],
                gene5primejunction=fusion["gene5prime_junction"],
                gene3prime=fusion["gene3prime"],
                gene3primejunction=fusion["gene3prime_junction"],
                db=db_human_hg38,
                pyensembl_data=data_human_hg38,
                protein_databases=["pfam"],
                noncanonical=False,
            )
            assert fusion.name in all_fusions, f"{fusion.name} not in list!"


if __name__ == "__main__":
    unittest.main()
