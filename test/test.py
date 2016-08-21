import agfusion

db = agfusion.AGFusionSQlite3DB('../bin/test.db','GRCh38')
#db.c.execute("SELECT * FROM GRCh38_transcript WHERE ensembl_gene_id==\"ENSG00000066468\"")
#t=db.c.fetchall()

fgfr2 = agfusion.Gene(
    query_gene="ENSG00000066468",
    junction=121485532,
    reference="GRCh38",
    db=db
)
fgfr2.predict_effect()

dnm3 = agfusion.Gene(
    query_gene="ENSG00000197959",
    junction=172253573,
    reference="GRCh38",
    db=db
)
dnm3.predict_effect()
