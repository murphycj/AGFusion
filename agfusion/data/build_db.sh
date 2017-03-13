
python ../../bin/agfusion_builddb --database agfusion.db --build homo_sapiens_core_84_38

python ../../bin/agfusion_builddb --database agfusion.db --build homo_sapiens_core_75_37

python ../../bin/agfusion_builddb --database agfusion.db --build mus_musculus_core_84_38

gzip agfusion.db
