import agfusion
import pyensembl
from Bio import SeqIO

def test_mouse_1():
    """
    test CDS and cDNA correct for junction that is on exon boundaries and
    produces an in-frame protein.
    """

    data = pyensembl.EnsemblRelease(84,'mouse')
    db = agfusion.AGFusionDB('../agfusion/data/agfusion.db')

    #test the dna and protein coding sequences are correct by comparing
    #with manually generally sequences

    fusion = agfusion.Fusion(
        gene5prime="ENSMUSG00000022770",
        gene5primejunction=31684294,
        gene3prime="ENSMUSG00000002413",
        gene3primejunction=39648486,
        db=db,
        pyensembl_data=data
    )

    fusion.save_transcript_cdna('DLG1-BRAF_mouse')
    fusion.save_transcript_cds('DLG1-BRAF_mouse')
    fusion.save_proteins('DLG1-BRAF_mouse')
    fusion.save_images('DLG1-BRAF_mouse')

    test_cdna = SeqIO.parse(open('Dlg1-Braf_cdna_manual.fa','r'),'fasta')
    test_cds = SeqIO.parse(open('Dlg1-Braf_cds_manual.fa','r'),'fasta')

    expected_transcript_combinations = [
        'ENSMUST00000100001-ENSMUST00000002487',
        'ENSMUST00000064477-ENSMUST00000002487',
        'ENSMUST00000115205-ENSMUST00000002487',
        'ENSMUST00000023454-ENSMUST00000002487',
        'ENSMUST00000115201-ENSMUST00000002487',
        'ENSMUST00000132176-ENSMUST00000002487'
    ]

    assert len(set(fusion.transcripts.keys()).intersection(set(expected_transcript_combinations)))==6, "Test 1: unexpected number protein coding transcripts."

    for seq in test_cdna:
        trans=fusion.transcripts[str(seq.id)]
        assert seq.seq==trans.cdna.seq, "cDNA is wrongly predicted: %s" % seq.id

    for seq in test_cds:
        trans=fusion.transcripts[str(seq.id)]
        assert seq.seq==trans.cds.seq, "cds is wrongly predicted: %s" % seq.id


def test_mouse_2():
    """
    Test CDS correct for junction within the exon (not on boundary) for two
    genes on reverse strand
    """

    data = pyensembl.EnsemblRelease(84,'mouse')
    db = agfusion.AGFusionDB('../agfusion/data/agfusion.db')

    fusion = agfusion.Fusion(
        gene5prime="ENSMUSG00000002413",
        gene5primejunction=39725110,
        gene3prime="ENSMUSG00000002413",
        gene3primejunction=39610402,
        db=db,
        pyensembl_data=data
    )

    cds = 'ATGGCGGCGCTGAGTGGCGGCGGTGGCAGCAGCAGCGGTGGCGGCGGCGGCGGTGGCGGCGGCGG' + \
          'TGGCGGTGGCGACGGCGGCGGCGGCGCCGAGCAGGGCCAGGCTCTGTTCAATGGCGACATGGAGC' + \
          'CGGAGGCCGGCGCTGGCGCCGCGGCCTCTTCGGCTGCGGACCCGGCCATTCCTGAAGAATTTGCAGCCTTCAAGTAG'

    assert str(fusion.transcripts['ENSMUST00000002487-ENSMUST00000002487'].cds.seq)==cds, "Test 2: CDS wrong"

def test_mouse_3():
    """
    Test CDS correct for junction within the exon (not on boundary) for one gene
    one the forward and one gene on the reverse strand
    """

    data = pyensembl.EnsemblRelease(84,'mouse')
    db = agfusion.AGFusionDB('../agfusion/data/agfusion.db')

    fusion = agfusion.Fusion(
        gene5prime="ENSMUSG00000022770",
        gene5primejunction=31664869,
        gene3prime="ENSMUSG00000002413",
        gene3primejunction=39610402,
        db=db,
        pyensembl_data=data
    )

    cds = 'ATGCCGGTCCGGAAGCAAGAATTTGCAGCCTTCAAGTAG'

    assert str(fusion.transcripts['ENSMUST00000064477-ENSMUST00000002487'].cds.seq)==cds, "Test 3: CDS wrong"

def test_mouse_4():
    """
    Test cDNA correctly produced for junctions being in UTRs
    """

    data = pyensembl.EnsemblRelease(84,'mouse')
    db = agfusion.AGFusionDB('../agfusion/data/agfusion.db')

    fusion = agfusion.Fusion(
        gene5prime="ENSMUSG00000022770",
        gene5primejunction=31664851,
        gene3prime="ENSMUSG00000022770",
        gene3primejunction=31873343,
        db=db,
        pyensembl_data=data
    )

    cdna = 'GGGGGTGCGGCCGCCGAAGGGGGAGCTCCTCCCCCGTCCCCTCACCCCCTCAGCTGAGCT' + \
          'CGGGGCGGGGCGGGGTACGTGGAGCGGGGCCGGGCGGGGAAGCTGCTCCGAGTCCGGCCG' + \
          'GAGCGCACCCGGGGCGCCCGCGTACGCCGCTCGCGGGAACTTTGCGGCGGAGCCGCAGGT' + \
          'GTGGAGGCCGCGGAGGGGGGTGCATGAGCGGCGCGGAGAGCGGCGGCTGTCCGGTCCGGC' + \
          'CCCTGCTGGAGTCGCCGCCGGGAGGAGACGAACGAGGAACCAG' + \
          'GTGTGTGCCGCCTTCCTGATTCTGGAGAAAA' + \
          'AAA'

    assert str(fusion.transcripts['ENSMUST00000064477-ENSMUST00000064477'].cdna.seq)==cdna, "Test 4: CDS wrong"

def test_mouse_5():
    """
    Test that AGFusion determines if the effect on each individual transcript
    """

    data = pyensembl.EnsemblRelease(84,'mouse')
    db = agfusion.AGFusionDB('../agfusion/data/agfusion.db')


    fusion = agfusion.Fusion(
        gene5prime="ENSMUSG00000022770",
        gene5primejunction=31665577,
        gene3prime="ENSMUSG00000002413",
        gene3primejunction=39651764,
        db=db,
        pyensembl_data=data
    )

    fusion = agfusion.Fusion(
        gene5prime="ENSMUSG00000022770",
        gene5primejunction=31664851,
        gene3prime="ENSMUSG00000002413",
        gene3primejunction=39610382,
        db=db,
        pyensembl_data=data
    )
    import pdb; pdb.set_trace()

    assert fusion.transcripts['ENSMUST00000064477-ENSMUST00000002487'].effect_5prime=='5UTR-end',"Test 5: Not found in 5'UTR-end"
    assert fusion.transcripts['ENSMUST00000064477-ENSMUST00000002487'].effect_3prime=='3UTR-start', "Test 5: Not found in at 3'UTR beginning"

    #All fusion isforms should be UTR-UTR

    fusion = agfusion.Fusion(
        gene5prime="ENSMUSG00000022770",
        gene5primejunction=31664851,
        gene3prime="ENSMUSG00000002413",
        gene3primejunction=39610382,
        db=db,
        pyensembl_data=data
    )

    assert fusion.transcripts['ENSMUST00000064477-ENSMUST00000002487'].effect_5prime=='5UTR-end',"Test 5: Not found in 5'UTR-end"
    assert fusion.transcripts['ENSMUST00000064477-ENSMUST00000002487'].effect_3prime=='3UTR-start', "Test 5: Not found in at 3'UTR beginning"

    fusion = agfusion.Fusion(
        gene5prime="ENSMUSG00000022770",
        gene5primejunction=31664850,
        gene3prime="ENSMUSG00000002413",
        gene3primejunction=39610382,
        db=db,
        pyensembl_data=data
    )
    assert fusion.transcripts['ENSMUST00000064477-ENSMUST00000002487'].effect_5prime=='5UTR',"Test 5: Not found in 5'UTR"
    assert fusion.transcripts['ENSMUST00000064477-ENSMUST00000002487'].effect_3prime=='3UTR-start', "Test 5: Not found in at 3'UTR beginning"


#test_mouse_1()
#test_mouse_2()
#test_mouse_3()
#test_mouse_4()
test_mouse_5()

print 'All tests passed!'
