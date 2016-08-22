import agfusion

db = agfusion.AGFusionDB('../../data/test.db')

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
