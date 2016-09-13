import agfusion
import pyensembl
from Bio import SeqIO

def test_mouse():
    data = pyensembl.EnsemblRelease(84,'mouse')

    db = agfusion.AGFusionDB('../data/agfusion.db',release=84,species='mouse')

    fgfr2 = agfusion.Gene(
        gene="Fgfr2",
        junction=130167703,
        db=db
    )

    dnm3 = agfusion.Gene(
        gene="Dnm3",
        junction=162076991,
        db=db
    )

    fusion = agfusion.model.Fusion(fgfr2,dnm3,db=db)
    fusion.save_transcript_cdna('fgfr2-dnm3_mouse')
    fusion.save_transcript_cds('fgfr2-dnm3_mouse')
    fusion.save_proteins('fgfr2-dnm3_mouse')
    #fusion.save_annotations('fgfr2-dnm3/domains.csv')
    #fusion.output_to_html()
    fusion.save_image('fgfr2-dnm3_mouse/fgfr2-dnm3')

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

    #fusion = agfusion.model.Fusion(fgfr2,dnm3,db=db,transcripts_5prime=['ENST00000346997'],transcripts_3prime=['ENST00000355305'])
    fusion = agfusion.model.Fusion(fgfr2,dnm3,db=db)
    fusion.save_transcript_cdna('fgfr2-dnm3')
    fusion.save_transcript_cds('fgfr2-dnm3')
    fusion.save_proteins('fgfr2-dnm3')
    #fusion.save_annotations('fgfr2-dnm3/domains.csv')
    fusion.save_image('fgfr2-dnm3')


    #assert cdna
    #cdna = SeqIO.parse(open('fgfr2-dnm3.cdna.fa','r'),'fasta').next()
    #trans=fusion.transcripts[u'ENST00000613048-ENST00000627582']
    #cdna.seq==trans.cdna.seq

test_mouse()
#test_human()
