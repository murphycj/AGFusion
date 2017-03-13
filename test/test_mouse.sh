../bin/agfusion \
  --gene5prime FGFR2 \
  --gene3prime DNM3 \
  --junction5prime 130167703 \
  --junction3prime 162019992 \
  --genome GRCm38 \
  --width 20 \
  --height 10 \
  --fontsize 8 \
  --out FGFR2-DNM3

../bin/agfusion \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF \
  --fontsize 12 \
  --height 3 \
  --width 8 \
  --dpi 90 \
  --colors Pkinase_Tyr:red L27_1:blue \
  --rename Pkinase_Tyr:Kinase L27_1:L27 \
  --WT

../bin/agfusion \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF-rescale \
  --colors "Serine-threonine/tyrosine-protein kinase catalytic domain;red" \
  --colors "L27_1;blue" \
  --rename "Serine-threonine/tyrosine-protein kinase catalytic domain;Kinase" \
  --rename "L27_1;L27" \
  --scale 2000
../bin/agfusion \
  --gene5prime FGFR2 \
  --gene3prime DNM3 \
  --junction5prime 130167703 \
  --junction3prime 162019992 \
  --genome GRCm38 \
  --out FGFR2-DNM3-rescale \
  --rename "Immunoglobulin I-set;I-set" \
  --rename "Dynamin GTPase effector;Dynamin" \
  --rename "Serine-threonine/tyrosine-protein kinase catalytic domain;Kinase" \
  --colors "Serine-threonine/tyrosine-protein kinase catalytic domain;red" \
  --scale 2000

../bin/agfusion \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF-noncanonical \
  --fontsize 12 \
  --height 3 \
  --width 8 \
  --dpi 90 \
  --colors Pkinase_Tyr:red L27_1:blue \
  --rename Pkinase_Tyr:Kinase L27_1:L27 \
  --dpi 100 \
  --noncanonical


../bin/agfusion \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF-reColorRename \
  --colors "Serine-threonine/tyrosine-protein kinase catalytic domain;red" \
  --colors "L27_1;blue" \
  --rename "Serine-threonine/tyrosine-protein kinase catalytic domain;Kinase" \
  --rename "L27_1;L27"
