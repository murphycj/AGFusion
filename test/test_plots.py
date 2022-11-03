"""
Tests for plotting.
"""
import unittest
from os.path import abspath, curdir, join
from pathlib import Path

import pyensembl

from agfusion import database, model

data = pyensembl.EnsemblRelease(84, "mouse")
db = database.AGFusionDB(abspath(join(curdir, "agfusion.mus_musculus.84.db")))
db.build = "mus_musculus_84"

data_human = pyensembl.EnsemblRelease(75, "human")
db_human = database.AGFusionDB(abspath(join(curdir, "agfusion.homo_sapiens.75.db")))
db_human.build = "homo_sapiens_75"


class TestSequencePredictionHuman(unittest.TestCase):
    """Test correctly predict human fusions"""

    def test_1(self):
        """
        test CDS and protein correct for junction that is on exon boundaries and
        produces an out-of-frame protein.
        """

        # test the dna and protein coding sequences are correct by comparing
        # with manually generally sequences

        fusion = model.Fusion(
            gene5prime="TMEM87B",
            gene5primejunction=112843681,
            gene3prime="MERTK",
            gene3primejunction=112722768,
            db=db_human,
            pyensembl_data=data_human,
            protein_databases=["pfam"],
            noncanonical=False,
        )

        fusion.save_images("DLG1-BRAF_mouse")

        assert Path(
            "DLG1-BRAF_mouse/ENST00000283206_ENST00000295408.png"
        ).exists(), "Could not save image."
