#
# Builds the agfusion databse for reference genomes GRCm38, GRCh38, and GRCh37
#

PERL5LIB=$PERL5LIB:~/Desktop/tools/perl-ensembl/BioPerl-1.6.1
PERL5LIB=$PERL5LIB:~/Desktop/tools/perl-ensembl/ensembl/modules/
PERL5LIB=$PERL5LIB:~/Desktop/tools/perl-ensembl/ensembl-compara/modules/
PERL5LIB=$PERL5LIB:~/Desktop/tools/perl-ensembl/ensembl-variation/modules/
PERL5LIB=$PERL5LIB:~/Desktop/tools/perl-ensembl/ensembl-funcgen/modules/
export PERL5LIB

perl ../../bin/builddb.pl -release 84 -genome GRCm38 -database agfusion.db -p 30

perl ../../bin/builddb.pl -release 84 -genome GRCh38 -database agfusion.db -p 30

unset PERL5LIB
PERL5LIB=$PERL5LIB:~/Desktop/tools/perl-ensembl/BioPerl-1.6.1
PERL5LIB=$PERL5LIB:~/Desktop/tools/ensembl-release-75/ensembl/modules/
PERL5LIB=$PERL5LIB:~/Desktop/tools/ensembl-release-75/ensembl-compara/modules/
PERL5LIB=$PERL5LIB:~/Desktop/tools/ensembl-release-75/ensembl-variation/modules/
PERL5LIB=$PERL5LIB:~/Desktop/tools/ensembl-release-75/ensembl-funcgen/modules/
export PERL5LIB

perl ../../bin/builddb.pl -release 75 -genome GRCh37 -database agfusion.db -p 30
