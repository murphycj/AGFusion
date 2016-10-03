import agfusion
import pyensembl
from Bio import SeqIO

def test_mouse():
    """
    This tests whether various gene fusions correspond to the expected output
    """

    data = pyensembl.EnsemblRelease(84,'mouse')
    db = agfusion.AGFusionDB('../data/agfusion.db')

    #test the dna and protein coding sequences are correct by comparing
    #with manually generally sequences

    dlg1 = agfusion.Gene(
        gene="ENSMUSG00000022770",
        junction=31684294,
        db=db,
        pyensembl_data=data
    )

    braf = agfusion.Gene(
        gene="ENSMUSG00000002413",
        junction=39648486,
        db=db,
        pyensembl_data=data
    )

    fusion = agfusion.model.Fusion(dlg1,braf,db=db,genome='GRCm38',middlestar=False)

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

    assert len(set(fusion.transcripts.keys()).intersection(set(expected_transcript_combinations)))==6, " unexpected number protein coding transcripts."

    for seq in test_cdna:
        trans=fusion.transcripts[str(seq.id)]
        assert seq.seq==trans.cdna.seq, "cDNA is wrongly predicted: %s" % seq.id

    for seq in test_cds:
        trans=fusion.transcripts[str(seq.id)]
        assert seq.seq==trans.cds.seq, "cds is wrongly predicted: %s" % seq.id

    #assert that the following fusions do not produce viable proteins
    #for various reasons

def test_human():

    data = pyensembl.EnsemblRelease(84,'human')
    db = agfusion.AGFusionDB('../data/agfusion.db')


    fgfr2 = agfusion.Gene(
        gene=data.gene_by_id("ENSG00000066468"),
        junction=121485532,
        reference="GRCh38",
        db=db
    )

    dnm3 = agfusion.Gene(
        gene=data.gene_by_id("ENSG00000197959"),
        junction=172407772,
        reference="GRCh38",
        db=db
    )

    fusion = agfusion.model.Fusion(fgfr2,dnm3,db=db)
    fusion.save_transcript_cdna('fgfr2-dnm3')
    fusion.save_transcript_cds('fgfr2-dnm3')
    fusion.save_proteins('fgfr2-dnm3')
    fusion.save_images('fgfr2-dnm3')


    #assert cdna
    #cdna = SeqIO.parse(open('fgfr2-dnm3.cdna.fa','r'),'fasta').next()
    #trans=fusion.transcripts[u'ENST00000613048-ENST00000627582']
    #cdna.seq==trans.cdna.seq

test_mouse()
#test_human()
