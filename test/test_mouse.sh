#../bin/agfusion \
#  --gene5prime FGFR2 \
#  --gene3prime DNM3 \
#  --junction5prime 130167703 \
#  --junction3prime 162019992 \
#  --genome GRCm38 \
#  --out FGFR2-DNM3 \
#  --colors Pkinase_Tyr:red \
#  --rename Pkinase_Tyr:Kinase \
#  --fontsize 12 \
#  --dpi 350 \
#  --height 2 \
#  --width 8


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
  --dpi 100
#  --scale 2000
