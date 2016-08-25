import agfusion
import pyensembl

data = pyensembl.EnsemblRelease(84,'human')

#db = agfusion.AGFusionDB('../../data/agfusion.db')

fgfr2 = agfusion.Gene(
    gene=data.gene_by_id("ENSG00000066468"),
    junction=121485532,
    reference="GRCh38",
    db=None
)
#fgfr2.predict_effect()

dnm3 = agfusion.Gene(
    gene=data.gene_by_id("ENSG00000197959"),
    junction=172253573,
    reference="GRCh38",
    db=None
)
#dnm3.predict_effect()

fusion = agfusion.model.Fusion(fgfr2,dnm3)
fusion.save_transcript_cdna('fgfr2-dnm3/transcript.cdna.fa')
fusion.save_transcript_cds('fgfr2-dnm3/transcript.cds.fa')
fusion.save_proteins('fgfr2-dnm3/protein.fa')
#fusion.save_image('fgfr2-dnm3')
