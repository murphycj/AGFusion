curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files//pfamA.txt.gz > 042117_pfamA.txt.gz

agfusion database --database agfusion.db --build homo_sapiens_core_84_38

agfusion database --database agfusion.db --build homo_sapiens_core_75_37

agfusion database --database agfusion.db --build mus_musculus_core_84_38

gzip agfusion.db
