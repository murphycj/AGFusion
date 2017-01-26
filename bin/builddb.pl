use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Getopt::ArgParse;
use DBD::SQLite;
use Parallel::ForkManager;
use Term::ProgressBar;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($INFO);

no warnings 'deprecated';

sub check_database_tables {
  my $dbh = shift;
  my $genome = shift;
  my $release = shift;

  #transcript table

  my $stmt = "CREATE TABLE if not exists "
      . "TRANSCRIPT_" . $genome . "_" . $release . " ("
      . "GENE_ID TEXT NOT NULL,"
      . "TRANSCRIPT_ID TEXT NOT NULL,"
      . "CANONICAL INTEGER)";

  my $rv = $dbh->do($stmt);

  if($rv < 0){
     print $DBI::errstr;
  }

  #protein features table

  my $stmt = "CREATE TABLE if not exists "
      . "PFEATURES_" . $genome . "_" . $release . " ("
      . "TRANSCRIPT_ID TEXT,"
      . "SOURCE TEXT,"
      . "ID TEXT,"
      . "DESCRIPTION TEXT,"
      . "SHORTDESCRIPTION TEXT,"
      . "START INTEGER,"
      . "END INTEGER)";

  my $rv = $dbh->do($stmt);

  if($rv < 0){
     print $DBI::errstr;
  }
}

sub insert_canonical {

  my $transcript = shift;
  my $dbh = shift;
  my $genome = shift;
  my $release = shift;

  my $canonical = $transcript->is_canonical();
  my $transcript_id = $transcript->stable_id();

  my $stmt = "INSERT INTO "
      . "TRANSCRIPT_" . $genome . "_" . $release
      . " (GENE_ID,TRANSCRIPT_ID,CANONICAL)"
      . " VALUES (" . $stable_id . "," . $transcript_id . "," . $canonical . ")";
  my $rv = $dbh->do($stmt) or die $DBI::errstr;
}

sub fetch_protein_features {
  #
  # fetch all the protein features associated with a transcript and return
  # an array
  #

  my $transcript = shift;
  my $transcript_id = $transcript->stable_id();
  my $dbh = shift;

  my $translation = $transcript->translation();

  my @features;

  if (defined $translation) {
    my $features = $translation->get_all_ProteinFeatures();

    while ( my $pfeature = shift @{$features} ) {
        my $source = $pfeature->analysis()->logic_name();
        my $id = $pfeature->display_id();
        my $desc = $pfeature->idesc();
        my $shortdesc = $pfeature->ilabel();
        my $start = $pfeature->start();
        my $end = $pfeature->end();

        push @features, [$transcript_id, $source, $id, $desc, $shortdesc, $start, $end];
    }
  }

  @features
}

$ap = Getopt::ArgParse->new_parser(
        prog        => 'MyProgramName',
        description => 'Build the SQLite3 database for a reference genomes by querying Ensembl perl API.'
);
$ap->add_arg(
  '-database',
  type => 'Scalar',
  required => 1,
  help => "Specify the location of the database."
);
$ap->add_arg(
  '-genome',
  type => 'Scalar',
  required => 1,
  help => "The Ensembl genome ID (GRCm38, GRCh38, or GRCh37)."
);
$ap->add_arg(
  '-release',
  type => 'Scalar',
  required => 1,
  help => "The Ensembl release number"
);
$ap->add_arg(
  '-p',
  type => 'Scalar',
  required => 1,
  help => "Number of processes."
);

my $args = $ap->parse_args();


#connect to database

INFO "Connecting to Ensembl servers...";

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

if ($args->genome=='GRCm38') {
  $gene_adaptor = $registry->get_adaptor('Mouse', 'Core', 'Gene');
} else {
  $gene_adaptor = $registry->get_adaptor('Human', 'Core', 'Gene');
}

my $dbh = DBI->connect("dbi:SQLite:dbname=" . $args->database,"","");

INFO "Checking database tables...";

check_database_tables($dbh,$args->genome,$args->release);

#loop over all genes

INFO "Fetching gene list from Ensembl...";

my $genes_ref = $gene_adaptor->fetch_all();
my @genes = @{$genes_ref};

my @gene_chunks;
push @gene_chunks, [ splice @genes, 0, 100 ] while @genes;

my @results;
my $counter = 0;

INFO "Fetching gene annotation (will take a while)...";

my $num_gene_chunks = scalar @gene_chunks;
my $progress = Term::ProgressBar->new($num_gene_chunks);

$pm = new Parallel::ForkManager($args->p);

$pm -> run_on_finish ( # called BEFORE the first call to start()
    sub {
      my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;

      $progress->update($_);

      push(@results, @$data_structure_reference);
    }
  );

while ($chunk = shift @gene_chunks) {

  my $pid = $pm->start and next;

  my @results_tmp;

  while ($gene = shift @{$chunk}) {

    #loop over the transcripts and insert into the table which are canonical
    #and then insert into the table the protein information

    my $transcripts = $gene->get_all_Transcripts();

    while ( my $transcript = shift @{$transcripts} ) {

        #insert_canonical($transcript,$dbh,$args->genome,$args->release);

        my @pfeatures = &fetch_protein_features($transcript,$dbh);
        push @results_tmp, @pfeatures;
    }
  }

  $pm->finish(0,\@results_tmp);

}

$pm->wait_all_children;

#add the data to the database

INFO "Adding fetched data to database...";

my $sth = $dbh->prepare(
  "INSERT INTO PFEATURES_" . $args->genome . "_" . $args->release
  . " (TRANSCRIPT_ID,SOURCE,ID,DESCRIPTION,SHORTDESCRIPTION,START,END) VALUES (?,?,?,?,?,?,?)"
);

foreach my $rec ( @results ) {
  $sth->execute( @$rec );
}
