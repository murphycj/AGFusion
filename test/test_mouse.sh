agfusion annotate \
  --gene5prime FGFR2 \
  --gene3prime DNM3 \
  --junction5prime 130167703 \
  --junction3prime 162019992 \
  -db agfusion.mus_musculus.84.db \
  --out FGFR2-DNM3

agfusion annotate \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  -db agfusion.mus_musculus.84.db \
  --out DLG1-BRAF

agfusion annotate \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --db agfusion.mus_musculus.84.db \
  --out DLG1-BRAF-reColorRename \
  --recolor "Pkinase_Tyr;red" --recolor "L27_1;blue" \
  --rename "Pkinase_Tyr;Kinase" --rename "L27_1;L27" \
  --WT

agfusion annotate \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --db agfusion.mus_musculus.84.db \
  --out DLG1-BRAF-rescale \
  --recolor "Pkinase_Tyr;red" --recolor "L27_1;blue" \
  --rename "Pkinase_Tyr;Kinase" --rename "L27_1;L27" \
  --scale 2000

agfusion annotate \
  --gene5prime FGFR2 \
  --gene3prime DNM3 \
  --junction5prime 130167703 \
  --junction3prime 162019992 \
  --db agfusion.mus_musculus.84.db \
  --out FGFR2-DNM3-rescale \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2000
