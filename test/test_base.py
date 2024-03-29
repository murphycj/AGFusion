"""Unit tests.
"""
import unittest
from os.path import abspath, curdir, join

import pyensembl
from Bio import Seq, SeqIO

from agfusion import database, model, parsers

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
            protein_databases=["pfam", "tmhmm"],
            noncanonical=False,
        )

        fusion.save_transcript_cdna("TMEM87B-MERTK-case0")
        fusion.save_transcript_cds("TMEM87B-MERTK-case0")
        fusion.save_proteins("TMEM87B-MERTK-case0")
        # fusion.save_images('DLG1-BRAF_mouse')

        test_cds = open("./data/test-human-case-0.txt", "r").read().strip()
        test_protein = Seq.Seq(test_cds).translate()
        test_protein = test_protein[0 : test_protein.find("*")]

        trans = fusion.transcripts["ENST00000283206_ENST00000295408"]

        assert test_cds == trans.cds.seq, "cds is wrongly predicted for human fusion (case 0)"
        assert (
            test_protein == trans.protein.seq
        ), "protein is wrongly predicted for human fusion (case 0)"

    def test_2(self):
        """
        test the dna and protein coding sequences are correct by comparing
        with manually generally sequences
        """

        fusion = model.Fusion(
            gene5prime="TMEM87B",
            gene5primejunction=112843681,
            gene3prime="MERTK",
            gene3primejunction=112722769,
            db=db_human,
            pyensembl_data=data_human,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=False,
        )

        fusion.save_transcript_cdna("TMEM87B-MERTK-case2")
        fusion.save_transcript_cds("TMEM87B-MERTK-case2")
        fusion.save_proteins("TMEM87B-MERTK-case2")
        # fusion.save_images('DLG1-BRAF_mouse')

        test_cds = open("./data/test-human-case-2.txt", "r").read().strip()
        test_protein = Seq.Seq(test_cds).translate()
        test_protein = test_protein[0 : test_protein.find("*")]

        trans = fusion.transcripts["ENST00000283206_ENST00000295408"]

        assert test_cds == trans.cds.seq, "cds is wrongly predicted for human fusion (case 2)"
        assert (
            test_protein == trans.protein.seq
        ), "protein is wrongly predicted for human fusion (case 2)"

    def test_3(self):
        """
        test the dna and protein coding sequences are correct by comparing
        with manually generally sequences
        """

        fusion = model.Fusion(
            gene5prime="TMEM87B",
            gene5primejunction=112843681,
            gene3prime="MERTK",
            gene3primejunction=112722771,
            db=db_human,
            pyensembl_data=data_human,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=False,
        )

        fusion.save_transcript_cdna("TMEM87B-MERTK-case3")
        fusion.save_transcript_cds("TMEM87B-MERTK-case3")
        fusion.save_proteins("TMEM87B-MERTK-case3")
        # fusion.save_images('DLG1-BRAF_mouse')

        test_cds = open("./data/test-human-case-3.txt", "r").read().strip()
        test_protein = Seq.Seq(test_cds).translate()
        test_protein = test_protein[0 : test_protein.find("*")]

        trans = fusion.transcripts["ENST00000283206_ENST00000295408"]

        assert test_cds == trans.cds.seq, "cds is wrongly predicted for human fusion (case 3)"
        assert (
            test_protein == trans.protein.seq
        ), "protein is wrongly predicted for human fusion (case 3)"


class TestSequencePrediction(unittest.TestCase):
    """Test correctly predict sequence content."""

    def test_1(self):
        """
        test CDS and cDNA correct for junction that is on exon boundaries and
        produces an in-frame protein.
        """

        # test the dna and protein coding sequences are correct by comparing
        # with manually generally sequences

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000022770",
            gene5primejunction=31684294,
            gene3prime="ENSMUSG00000002413",
            gene3primejunction=39648486,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )

        fusion.save_transcript_cdna("DLG1-BRAF_mouse")
        fusion.save_transcript_cds("DLG1-BRAF_mouse")
        fusion.save_proteins("DLG1-BRAF_mouse")
        # fusion.save_images('DLG1-BRAF_mouse')

        test_cdna = SeqIO.parse(open("./data/Dlg1-Braf_cdna_manual.fa", "r"), "fasta")
        test_cds = SeqIO.parse(open("./data/Dlg1-Braf_cds_manual.fa", "r"), "fasta")

        expected_transcript_combinations = [
            "ENSMUST00000100001_ENSMUST00000002487",
            "ENSMUST00000064477_ENSMUST00000002487",
            "ENSMUST00000115205_ENSMUST00000002487",
            "ENSMUST00000023454_ENSMUST00000002487",
            "ENSMUST00000115201_ENSMUST00000002487",
            "ENSMUST00000132176_ENSMUST00000002487",
        ]

        assert (
            len(set(fusion.transcripts.keys()).intersection(set(expected_transcript_combinations)))
            == 6
        ), "Test 1: unexpected number protein coding transcripts."

        for seq in test_cdna:
            trans = fusion.transcripts[str(seq.id)]
            assert seq.seq == trans.cdna.seq, f"cDNA is wrongly predicted: {seq.id}"

        for seq in test_cds:
            trans = fusion.transcripts[str(seq.id)]
            assert seq.seq == trans.cds.seq, f"cds is wrongly predicted: {seq.id}"

    def test_2(self):
        """
        Test CDS correct for junction within the exon (not on boundary) for two
        genes on reverse strand
        """

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000002413",
            gene5primejunction=39725110,
            gene3prime="ENSMUSG00000002413",
            gene3primejunction=39610402,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )

        cds = (
            "ATGGCGGCGCTGAGTGGCGGCGGTGGCAGCAGCAGCGGTGGCGGCGGCGGCGGTGGCGGCGGCGG"
            + "TGGCGGTGGCGACGGCGGCGGCGGCGCCGAGCAGGGCCAGGCTCTGTTCAATGGCGACATGGAGC"
            + "CGGAGGCCGGCGCTGGCGCCGCGGCCTCTTCGGCTGCGGACCCGGCCATTCCTGAAGAATTTGCAGCCTTCAAGTAG"
        )

        assert (
            str(fusion.transcripts["ENSMUST00000002487_ENSMUST00000002487"].cds.seq) == cds
        ), "Test 2: CDS wrong"

    def test_3(self):
        """
        Test CDS correct for junction within the exon (not on boundary) for one gene
        one the forward and one gene on the reverse strand
        """

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000022770",
            gene5primejunction=31664869,
            gene3prime="ENSMUSG00000002413",
            gene3primejunction=39610402,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )

        cds = "ATGCCGGTCCGGAAGCAAGAATTTGCAGCCTTCAAGTAG"

        assert (
            str(fusion.transcripts["ENSMUST00000064477_ENSMUST00000002487"].cds.seq) == cds
        ), "Test 3: CDS wrong"

    def test_mouse_4(self):
        """
        Test cDNA correctly produced for junctions being in UTRs
        """

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000022770",
            gene5primejunction=31664851,
            gene3prime="ENSMUSG00000022770",
            gene3primejunction=31873343,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )

        cdna = (
            "GGGGGTGCGGCCGCCGAAGGGGGAGCTCCTCCCCCGTCCCCTCACCCCCTCAGCTGAGCT"
            + "CGGGGCGGGGCGGGGTACGTGGAGCGGGGCCGGGCGGGGAAGCTGCTCCGAGTCCGGCCG"
            + "GAGCGCACCCGGGGCGCCCGCGTACGCCGCTCGCGGGAACTTTGCGGCGGAGCCGCAGGT"
            + "GTGGAGGCCGCGGAGGGGGGTGCATGAGCGGCGCGGAGAGCGGCGGCTGTCCGGTCCGGC"
            + "CCCTGCTGGAGTCGCCGCCGGGAGGAGACGAACGAGGAACCAG"
            + "GTGTGTGCCGCCTTCCTGATTCTGGAGAAAA"
            + "AAAA"
        )

        assert (
            str(fusion.transcripts["ENSMUST00000064477_ENSMUST00000064477"].cdna.seq) == cdna
        ), "Test 4: cDNA wrong"

    def test_mouse_5(self):
        """Test transcript constructed correctly (mouse)"""

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000022770",
            gene5primejunction=31664850,
            gene3prime="ENSMUSG00000002413",
            gene3primejunction=39603240,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )

        cdna = (
            "GGGGGTGCGGCCGCCGAAGGGGGAGCTCCTCCCCCGTCCCCTCACCCCCTCAGCTGAGCT"
            + "CGGGGCGGGGCGGGGTACGTGGAGCGGGGCCGGGCGGGGAAGCTGCTCCGAGTCCGGCCG"
            + "GAGCGCACCCGGGGCGCCCGCGTACGCCGCTCGCGGGAACTTTGCGGCGGAGCCGCAGGT"
            + "GTGGAGGCCGCGGAGGGGGGTGCATGAGCGGCGCGGAGAGCGGCGGCTGTCCGGTCCGGC"
            + "CCCTGCTGGAGTCGCCGCCGGGAGGAGACGAACGAGGAACCAG"
            + "GTGTGTGCCGCCTTCCTGATTCTGGAGAAA"
            + "GAAA"
        )

        assert (
            str(fusion.transcripts["ENSMUST00000064477_ENSMUST00000002487"].cdna.seq) == cdna
        ), "Test 10: incorrect cDNA"


class TestEffectPrediciton(unittest.TestCase):
    """Test correctly predict junction location annotation and fusion effect."""

    def test_mouse_1(self):
        """
        Test that AGFusion determines if the effect on each individual transcript
        """

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000022770",
            gene5primejunction=31664852,
            gene3prime="ENSMUSG00000002413",
            gene3primejunction=39651764,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )

        assert (
            fusion.transcripts["ENSMUST00000064477_ENSMUST00000002487"].effect_5prime
            == "CDS (start)"
        ), "Test 5: not CDS start"
        assert (
            fusion.transcripts["ENSMUST00000064477_ENSMUST00000002487"].effect_3prime == "CDS"
        ), "Test 5: not CDS"

    def test_mouse_2(self):
        """Test junction location annotated"""

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000022770",
            gene5primejunction=31664851,
            gene3prime="ENSMUSG00000002413",
            gene3primejunction=39610381,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )

        assert (
            fusion.transcripts["ENSMUST00000064477_ENSMUST00000002487"].effect_5prime
            == "5UTR (end)"
        ), "Test 6: Not found in 5'UTR-end"
        assert (
            fusion.transcripts["ENSMUST00000064477_ENSMUST00000002487"].effect_3prime
            == "3UTR (start)"
        ), "Test 6: Not found in at 3'UTR beginning"

    def test_mouse_3(self):
        """Test junction location annotated"""

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000022770",
            gene5primejunction=31664850,
            gene3prime="ENSMUSG00000002413",
            gene3primejunction=39610381,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )
        assert (
            fusion.transcripts["ENSMUST00000064477_ENSMUST00000002487"].effect_5prime == "5UTR"
        ), "Test 7: Not found in 5'UTR"
        assert (
            fusion.transcripts["ENSMUST00000064477_ENSMUST00000002487"].effect_3prime
            == "3UTR (start)"
        ), "Test 7: Not found in at 3'UTR beginning"

    def test_mouse_4(self):
        """Test junction location annotated"""

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000022770",
            gene5primejunction=31664851,
            gene3prime="ENSMUSG00000022770",
            gene3primejunction=31871782,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )

        assert (
            fusion.transcripts["ENSMUST00000064477_ENSMUST00000064477"].effect_5prime
            == "5UTR (end)"
        ), "Test 8: Not found in 5'UTR-end"
        assert (
            fusion.transcripts["ENSMUST00000064477_ENSMUST00000064477"].effect_3prime
            == "3UTR (start)"
        ), "Test 8: Not found in at 3'UTR beginning"

    def test_mouse_5(self):
        """Test junction location annotated"""

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000002413",
            gene5primejunction=39725296,
            gene3prime="ENSMUSG00000002413",
            gene3primejunction=39610381,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )

        assert (
            fusion.transcripts["ENSMUST00000002487_ENSMUST00000002487"].effect_5prime
            == "5UTR (end)"
        ), "Test 9: Not found in 5'UTR-end"
        assert (
            fusion.transcripts["ENSMUST00000002487_ENSMUST00000002487"].effect_3prime
            == "3UTR (start)"
        ), "Test 9: Not found in at 3'UTR beginning"

    def test_mouse_6(self):
        """Test junction location annotated"""

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000022770",
            gene5primejunction=31743271,
            gene3prime="ENSMUSG00000002413",
            gene3primejunction=39665003,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )
        transcript = fusion.transcripts["ENSMUST00000023454_ENSMUST00000002487"]

        assert (
            transcript.effect_5prime == "intron (cds)"
        ), f"Test 11: incorrect 5' effect: {transcript.effect_5prime}"
        assert (
            transcript.effect_3prime == "intron (cds)"
        ), f"Test 11: incorrect 3' effect: {transcript.effect_3prime}"

    def test_mouse_7(self):
        """Test junction location annotated"""

        fusion = model.Fusion(
            gene5prime="ENSMUSG00000022770",
            gene5primejunction=31664820,
            gene3prime="ENSMUSG00000002413",
            gene3primejunction=39610405,
            db=db,
            pyensembl_data=data,
            protein_databases=["pfam", "tmhmm"],
            noncanonical=True,
        )
        transcript = fusion.transcripts["ENSMUST00000023454_ENSMUST00000002487"]

        assert (
            transcript.effect_5prime == "intron (before cds)"
        ), f"Test 12: incorrect 5' effect: {transcript.effect_5prime}"
        assert (
            transcript.effect_3prime == "intron (cds)"
        ), f"Test 12: incorrect 3' effect: {transcript.effect_3prime}"


class TestFusionCatcher(unittest.TestCase):
    """Test parse FusionCatcher parse."""

    def test_1(self):
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
            "./data/FusionsFindingAlgorithms/FusionCatcher/final-list_candidate-fusion-genes.txt",
            db.logger,
        ):
            fusion = model.Fusion(
                gene5prime=fusion["gene5prime"],
                gene5primejunction=fusion["gene5prime_junction"],
                gene3prime=fusion["gene3prime"],
                gene3primejunction=fusion["gene3prime_junction"],
                db=db,
                pyensembl_data=data,
                protein_databases=["pfam"],
                noncanonical=False,
            )
            assert fusion.name in all_fusions, f"{fusion.name} not in list!"


class TestSTARFusion(unittest.TestCase):
    """Test parse STAR Fusion"""

    def test_1(self):
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
            "./data/FusionsFindingAlgorithms/FusionCatcher/final-list_candidate-fusion-genes.txt",
            db.logger,
        ):
            fusion = model.Fusion(
                gene5prime=fusion["gene5prime"],
                gene5primejunction=fusion["gene5prime_junction"],
                gene3prime=fusion["gene3prime"],
                gene3primejunction=fusion["gene3prime_junction"],
                db=db,
                pyensembl_data=data,
                protein_databases=["pfam"],
                noncanonical=False,
            )
            assert fusion.name in all_fusions, f"{fusion.name} not in list!"


if __name__ == "__main__":
    unittest.main()
