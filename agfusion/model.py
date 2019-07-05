import itertools
import os
import re
import sys

from agfusion import utils, exceptions, plot
import pandas
from Bio import Seq, SeqIO, SeqRecord, SeqUtils
from Bio.Alphabet import generic_dna, generic_protein

# matplotlib.rcParams['interactive'] = False

from agfusion.utils import STANDARD_CHROMOSOMES, MIN_DOMAIN_LENGTH

class _Gene():
    """
    Stores the necessary information to specify the architecture of either
    wild-type genes or fusion gene.
    """

    def __init__(self, genes=None, junction=0, pyensembl_data=None,
                 genome='', gene5prime=False, db=None, noncanonical=False):
        """
        genes : str or list
            Provide one gene (str) or list of genes (list). In the case of a
            list of genes, it will seach each gene until it finds a hit.

        junction : int

        pyensembl_data : None

        genome : str

        gene5prime : bool

        db : None

        noncanonical : bool
        """

        if type(genes)==str:
            genes = [genes]
        elif type(genes)!=list:
            db.logger.error(
                'parameter \'gene\' to _Gene is not string or list: {}'
                .format(genes)
            )
            exit()
        if type(junction)!=int:
            db.logger.error(
                'parameter \'junction\' to _Gene is not int: {}'
                .format(junction)
            )
            exit()

        self.gene_found = False
        self.provided_transcript = False # indicates user if user provided transcript

        self.domains = []
        self.transcripts = {} # stores transcripts and their DB keys
        self.gene = None
        self.pyensembl_data = pyensembl_data
        self.genome = genome
        self.gene5prime = gene5prime
        self.db = db
        self.noncanonical = noncanonical

        # find the appropriate Ensembl gene ID

        for gene in genes:

            if re.findall('^NM_',gene) or re.findall('^NR_',gene):
                self._search_by_refseq(gene)
            if gene.isdigit() and not self.gene_found:
                self._search_as_entrez(gene)
            if re.findall('(^ENS.*G)', gene.upper()) and not self.gene_found:
                self._search_as_ensembl_id(gene)
            if re.findall('(^ENS.*T)', gene.upper()) and not self.gene_found:
                self._search_as_ensembl_transcript_id(gene)

            if not self.gene_found:

                # else check if it is a gene symbol

                self._search_by_symbol(gene)

                if not self.gene_found:
                    gene = gene.capitalize()
                    self._search_by_symbol(gene)

                if not self.gene_found:
                    gene = gene.upper()
                    self._search_by_symbol(gene)

            if self.gene_found:
                break

        # if gene has not been identified yet

        if not self.gene_found:
            raise exceptions.GeneIDException(','.join(genes))

        # else continue with processing

        self.junction = junction

        if not self.gene.contains(self.gene.contig,self.junction,self.junction):
            raise exceptions.JunctionException(self.gene.name, self.junction)

        # fetch the entrez gene id and canonical transcript id

        sqlite3_command = "SELECT * FROM " + self.db.build + " WHERE stable_id==\"" + self.gene.gene_id + "\""
        self.db.logger.debug('SQLite - ' + sqlite3_command)
        self.db.sqlite3_cursor.execute(
            sqlite3_command
        )
        tmp = self.db.sqlite3_cursor.fetchall()

        self.gene_id = tmp[0][0]
        self.entrez_id = tmp[0][2]
        self.canonical_transcript_id = tmp[0][4]

        # get the transcripts that will be processed and annotated

        if not noncanonical and not self.provided_transcript:

            # if only want the canonical and did not specify a certain transcript

            sqlite3_command = "SELECT * FROM " + self.db.build + "_transcript WHERE transcript_id==\"" + self.canonical_transcript_id + "\""
            self.db.logger.debug('SQLite - ' + sqlite3_command)
            self.db.sqlite3_cursor.execute(
                sqlite3_command
            )
            tmp = db.sqlite3_cursor.fetchall()
            self.transcripts[tmp[0][2]] = tmp[0][0]
        elif not noncanonical and self.provided_transcript:
            pass
        else:
            if self.provided_transcript and noncanonical:
                self.db.logger.warn("You provided a transcript ID as well as specified --noncanonical flag. Will process all the gene's transcripts.")

            # fetch transcript ids

            sqlite3_command = "SELECT * FROM " + self.db.build + "_transcript WHERE gene_id==\"" + self.gene_id + "\""
            self.db.logger.debug('SQLite - ' + sqlite3_command)
            self.db.sqlite3_cursor.execute(
                sqlite3_command
            )
            tmp = self.db.sqlite3_cursor.fetchall()

            for transcript in tmp:
                self.transcripts[transcript[2]] = transcript[0]

    def _search_as_ensembl_transcript_id(self,gene):
        # if it is ensembl transcript id

        if gene in self.pyensembl_data.transcript_ids():
            transcript = self.pyensembl_data.transcript_by_id(gene.upper())
            self.gene = transcript.gene

            self.db.logger.debug('Found Ensembl transcript entry for %s: %s' % (gene,self.gene.id))
            self.gene_found = True
            self.provided_transcript = True
        else:
            self.db.logger.debug('Found no Ensembl transcript entry for %s' % gene)

        sqlite3_command = "SELECT * FROM " + self.db.build + "_transcript WHERE transcript_stable_id==\"" + gene + "\""
        self.db.logger.debug('SQLite - ' + sqlite3_command)
        self.db.sqlite3_cursor.execute(
            sqlite3_command
        )
        tmp = self.db.sqlite3_cursor.fetchall()
        self.transcripts[tmp[0][2]] = tmp[0][0]

    def _search_as_ensembl_id(self,gene):
        # if it is ensembl gene id

        if gene in self.pyensembl_data.gene_ids():
            self.gene = self.pyensembl_data.gene_by_id(gene.upper())
            self.db.logger.debug('Found Enrez gene ID entry for %s: %s' % (gene,self.gene.id))
            self.gene_found = True
        else:
            self.db.logger.debug('Cannot find Ensembl gene id %s in database!' % gene)

    def _search_as_entrez(self,gene):
        # if it is an entrez gene ID

        sqlite3_command = "SELECT * FROM " + self.db.build + " WHERE entrez_id==\"" + gene + "\""
        self.db.logger.debug('SQLite - ' + sqlite3_command)
        self.db.sqlite3_cursor.execute(
            sqlite3_command
        )
        tmp = self.db.sqlite3_cursor.fetchall()

        if len(tmp)==1:
            self.gene = self.pyensembl_data.gene_by_id(tmp[0][1])
            self.db.logger.debug('Found Entrez gene ID entry for %s: %s' % (gene,self.gene.id))
            self.gene_found = True
        elif len(tmp)>1:

            self.db.logger.error('Found too many Entrez gene ID entries for %s!' % gene)
            for i in tmp:
                self.db.logger.error(i)
            sys.exit()
        else:
            self.db.logger.debug('Found no Entrez gene ID entry for %s' % gene)

    def _search_by_refseq(self,gene):


        # if it is RefSeq ID

        sqlite3_command = "SELECT * FROM " + self.db.build + "_refseq WHERE refseq_id==\"" + gene + "\""
        self.db.logger.debug('SQLite - ' + sqlite3_command)
        self.db.sqlite3_cursor.execute(
            sqlite3_command
        )
        tmp = self.db.sqlite3_cursor.fetchall()

        if len(tmp)==1:
            self.transcripts[tmp[0][1]] = tmp[0][0]
            self.db.logger.debug('Found RefSeq entry for %s' % gene)
            self.gene_found = True
            self.provided_transcript = True

            #fetch the gene

            transcript = self.pyensembl_data.transcript_by_id(tmp[0][1])
            self.gene = transcript.gene

        elif len(tmp)>1:
            self.db.logger.error('Found too many RefSeq entries for %s!' % gene)
            for i in tmp:
                self.db.logger.error(i)
            sys.exit()
        else:
            self.db.logger.debug('Found no RefSeq entry for %s' % gene)

    def _search_by_symbol(self,gene):
        if gene in self.pyensembl_data.gene_names():
            temp = self.pyensembl_data.genes_by_name(gene)

            if len(temp) > 1:

                # if too many ensembl gene IDs returned
                # use the one located on chromosomes

                for tmp_gene in temp:
                    if tmp_gene.contig in STANDARD_CHROMOSOMES:
                        self.gene = tmp_gene
                        self.gene_found = True

                if self.gene is None:
                    raise exceptions.TooManyGenesException(
                        gene,
                        [x.id for x in temp],
                        self.genome
                    )
            else:
                self.gene = temp[0]
                self.gene_found = True
                self.db.logger.debug('Found gene symbol entry for %s: %s' % (gene,self.gene.id))

class Fusion():
    """
    Generates the information needed for the gene fusion
    """

    def __init__(
            self,
            gene5prime=None, gene5primejunction=0,
            gene3prime=None, gene3primejunction=0,
            db=None, pyensembl_data=None, protein_databases=None,
            noncanonical=False,
            transcripts_5prime=None, transcripts_3prime=None):
        """
        gene5prime : str

        gene3prime : str

        gene5primejunction : int

        gene3primejunction : int

        db : None

        pyensembl_data : None

        protein_databases : None

        noncanonical : bool

        transcripts_5prime : str

        transcripts_3prime : str
        """

        self.db = db
        self.pyensembl_data = pyensembl_data

        # get the reference genom

        self.gene5prime = _Gene(
            genes=gene5prime,
            junction=gene5primejunction,
            pyensembl_data=pyensembl_data,
            gene5prime=True,
            db=db,
            noncanonical=noncanonical
        )

        self.gene3prime = _Gene(
            genes=gene3prime,
            junction=gene3primejunction,
            pyensembl_data=pyensembl_data,
            gene5prime=False,
            db=db,
            noncanonical=noncanonical
        )

        self.name = self.gene5prime.gene.name + '-' + self.gene3prime.gene.name

        # fetch which are canonical transcripts

        gene5prime_canonical = ''
        gene3prime_canonical = ''

        # construct all the fusion transcript combinations

        self.transcripts = {}

        for combo in list(itertools.product(
                        list(self.gene5prime.transcripts.keys()),
                        list(self.gene3prime.transcripts.keys()))):
            transcript1 = pyensembl_data.transcript_by_id(combo[0])
            transcript2 = pyensembl_data.transcript_by_id(combo[1])

            if transcripts_5prime is not None \
                    and transcript1.id not in transcripts_5prime:

                continue

            if transcripts_3prime is not None \
                    and transcript2.id not in transcripts_3prime:

                continue

            # skip if the junction is outside the range of either transcript

            name = transcript1.id + '-' + transcript2.id

            if not transcript1.contains(transcript1.contig,self.gene5prime.junction,self.gene5prime.junction):
                self.transcripts[name] = FusionTranscript(
                    transcript1=transcript1,
                    transcript2=transcript2,
                    gene5prime=self.gene5prime,
                    gene3prime=self.gene3prime,
                    db=db,
                    protein_databases=protein_databases,
                    predict_effect=False
                )
                self.transcripts[name].effect = 'Outside transcript boundry'
                self.transcripts[name].has_coding_potential = False

            elif not transcript2.contains(transcript2.contig,self.gene3prime.junction,self.gene3prime.junction):
                self.transcripts[name] = FusionTranscript(
                    transcript1=transcript1,
                    transcript2=transcript2,
                    gene5prime=self.gene5prime,
                    gene3prime=self.gene3prime,
                    db=db,
                    protein_databases=protein_databases,
                    predict_effect=False
                )
                self.transcripts[name].effect='Outside transcript boundry'
                self.transcripts[name].has_coding_potential=False

            else:
                self.transcripts[name] = FusionTranscript(
                    transcript1=transcript1,
                    transcript2=transcript2,
                    gene5prime=self.gene5prime,
                    gene3prime=self.gene3prime,
                    db=db,
                    protein_databases=protein_databases,
                )

    def output_to_html(
            self, fontsize=12, dpi=90, colors={}, rename={},
            width=8, height=2, scale=0, mpld3=None):
        """
        Write figures into html format for plotting in web tool rather
        than writing figures to file

        Needs the mpld3 package to be used
        """

        dict_of_plots = list()
        plot_key = dict()

        for name, transcript in list(self.transcripts.items()):

            if not transcript.has_coding_potential:
                continue

            fig = plt.figure(figsize=(width, height), dpi=dpi, frameon=False)

            self._draw(
                fig=fig,
                name=name,
                transcript=transcript,
                scale=scale,
                fontsize=fontsize,
                colors=[],
                rename=[]
            )

            single_chart = dict()
            single_chart['id'] = name
            single_chart['json'] = json.dumps(mpld3.fig_to_dict(fig))
            dict_of_plots.append(single_chart)

            plot_key[name] = transcript.name

            plt.close(fig)
            plt.clf()

        return dict_of_plots, plot_key

    def save_images(
            self, out_dir='', file_type='png', fontsize=12, dpi=100,
            colors={}, rename={}, width=8, height=2, scale=0,
            no_domain_labels=True, plot_WT=False,exclude=None):
        """
        Save images of all fusion isoforms and write to file
        Also plot the wild type proteins if desired
        """

        assert file_type in ['png', 'pdf'], 'provided wrong file type'

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        if plot_WT:
            gene5prime_WT = os.path.join(
                out_dir,
                self.gene5prime.gene.gene_name
            )
            gene3prime_WT = os.path.join(
                out_dir,
                self.gene3prime.gene.gene_name
            )
            if not os.path.exists(gene5prime_WT):
                os.mkdir(gene5prime_WT)

            if not os.path.exists(gene3prime_WT):
                os.mkdir(gene3prime_WT)

        for name, transcript in list(self.transcripts.items()):

            if not transcript.has_coding_potential:
                continue

            filename = os.path.join(out_dir, name + '.' + file_type)

            pplot = plot.PlotFusionProtein(
                filename=filename,
                width=width,
                height=height,
                dpi=dpi,
                scale=scale,
                fontsize=fontsize,
                colors=colors,
                rename=rename,
                no_domain_labels=no_domain_labels,
                transcript=transcript,
                exclude=exclude
            )
            pplot.draw()
            pplot.save()

            filename = os.path.join(
                out_dir,
                name + '.exon.' + file_type
            )

            pplot = plot.PlotFusionExons(
                transcript=transcript,
                filename=filename,
                width=width,
                height=height,
                dpi=dpi,
                scale=scale,
                fontsize=fontsize
            )
            pplot.draw()
            pplot.save()

            if plot_WT:

                # plot exons

                filename = os.path.join(
                    gene5prime_WT,
                    transcript.transcript1.id + '.exon.' + file_type
                )

                pplot = plot.PlotWTExons(
                    ensembl_transcript=transcript.transcript1,
                    filename=filename,
                    width=width,
                    height=height,
                    dpi=dpi,
                    scale=scale,
                    fontsize=fontsize
                )
                pplot.draw()
                pplot.save()

                filename = os.path.join(
                    gene3prime_WT,
                    transcript.transcript2.id + '.exon.' + file_type
                )

                pplot = plot.PlotWTExons(
                    ensembl_transcript=transcript.transcript2,
                    filename=filename,
                    width=width,
                    height=height,
                    dpi=dpi,
                    scale=scale,
                    fontsize=fontsize
                )
                pplot.draw()
                pplot.save()

                # plot proteins

                filename = os.path.join(
                    gene5prime_WT,
                    transcript.transcript1.id + '.' + file_type
                )

                pplot = plot.PlotWTProtein(
                    ensembl_transcript=transcript.transcript1,
                    filename=filename,
                    width=width,
                    height=height,
                    dpi=dpi,
                    scale=scale,
                    fontsize=fontsize,
                    colors=colors,
                    rename=rename,
                    no_domain_labels=no_domain_labels,
                    transcript=transcript,
                    exclude=exclude
                )
                pplot.draw()
                pplot.save()

                filename = os.path.join(
                    gene3prime_WT,
                    transcript.transcript2.id + '.' + file_type
                )

                pplot = plot.PlotWTProtein(
                    ensembl_transcript=transcript.transcript2,
                    filename=filename,
                    width=width,
                    height=height,
                    dpi=dpi,
                    scale=scale,
                    fontsize=fontsize,
                    colors=colors,
                    rename=rename,
                    no_domain_labels=no_domain_labels,
                    transcript=transcript,
                    exclude=exclude
                )
                pplot.draw()
                pplot.save()

    def save_transcript_cdna(self, out_dir='.', middlestar=False):
        """
        Save the cDNA sequences for all fusion isoforms to a fasta file
        """

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        fout = open(
            os.path.join(
                out_dir,
                self.name + '_cdna.fa'
            ),
            'w'
        )

        for name, transcript in list(self.transcripts.items()):

            if transcript.cdna is not None:

                if middlestar:
                    temp = str(transcript.cdna.seq)
                    temp = temp[:transcript.transcript_cdna_junction_5prime] + '*' + temp[transcript.transcript_cdna_junction_5prime:]
                    transcript.cdna.seq = Seq.Seq(temp,generic_dna)

                SeqIO.write(transcript.cdna,fout,"fasta")
            else:
                cdna = SeqRecord.SeqRecord(
                    Seq.Seq("",generic_dna),
                    id=transcript.name,
                    name=transcript.name,
                    description="No cDNA, fusion junction outside transcript(s) boundary"
                )
                SeqIO.write(cdna,fout,"fasta")

        fout.close()

    def save_transcript_cds(self, out_dir='.', middlestar=False):
        """
        Save the CDS sequences for all fusion isoforms to a fasta file
        """

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        #check if any transcripts have coding potential

        n=0
        for name, transcript in list(self.transcripts.items()):

            if transcript.cds is not None:
                n+=1

        if n == 0:
            self.db.logger.debug('The %s fusion does not produce any protein coding transcripts. No cds.fa file will be written' % self.name)
            return

        fout = open(
            os.path.join(
                out_dir,
                self.name + '_cds.fa'
            ),
            'w'
        )

        for name, transcript in list(self.transcripts.items()):

            if transcript.cds is not None:

                if middlestar:
                    temp = str(transcript.cds.seq)
                    temp = temp[:transcript.transcript_cds_junction_5prime] + '*' + temp[transcript.transcript_cds_junction_5prime:]
                    transcript.cds.seq = Seq.Seq(temp,generic_dna)

                SeqIO.write(transcript.cds,fout,"fasta")

        fout.close()

    def save_proteins(self, out_dir='.', middlestar=False):
        """
        Save the protein sequences for all fusion isoforms to a fasta file
        """

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        # check if any transcripts have coding potential

        n = 0
        for name, transcript in list(self.transcripts.items()):

            if transcript.cds is not None:
                n += 1

        if n == 0:
            self.db.logger.debug('The %s fusion does not produce any protein coding transcripts. No proteins.fa file will be written' % self.name)
            return

        fout = open(
            os.path.join(
                out_dir,
                self.name + '_protein.fa'
            ),
            'w'
        )

        for name, transcript in list(self.transcripts.items()):

            if transcript.cds is not None:

                if middlestar:
                    temp = str(transcript.protein.seq)
                    temp = temp[:transcript.transcript_protein_junction_5prime] + '*' + temp[transcript.transcript_protein_junction_5prime:]
                    transcript.protein.seq = Seq.Seq(temp,generic_protein)

                SeqIO.write(transcript.protein,fout,"fasta")

        fout.close()

    def save_tables(self,out_dir='.',annotation='pfam'):

        #write table of fusion isoforms

        fout = open(
            os.path.join(
                out_dir,
                self.name + '.fusion_transcripts.txt'
            ),
            'w'
        )
        fout.write(
            ','.join(['{}']*11).format(
                "5\'_gene",
                "3\'_gene",
                "5\'_transcript",
                "3\'_transcript",
                "5\'_strand",
                "3\'_strand",
                "5\'_transcript_biotype",
                "3\'_transcript_biotype",
                "Fusion_effect",
                "Protein_length",
                "Protein_weight_(kD)"
            ) + '\n')

        for name, transcript in list(self.transcripts.items()):

            if transcript.protein_length is None:
                protein_length = "NA"
            else:
                protein_length = transcript.protein_length

            if transcript.molecular_weight is None:
                molecular_weight = "NA"
            else:
                molecular_weight = transcript.molecular_weight

            fout.write(
                ','.join(['{}']*11).format(
                    transcript.gene5prime.gene.gene_name,
                    transcript.gene3prime.gene.gene_name,
                    transcript.transcript1.id,
                    transcript.transcript2.id,
                    transcript.transcript1.strand,
                    transcript.transcript2.strand,
                    transcript.transcript1.biotype,
                    transcript.transcript2.biotype,
                    transcript.effect,
                    protein_length,
                    molecular_weight))
        fout.close()

        # write table containing protein domains for each fusion isoform

        fout = open(
            os.path.join(
                out_dir,
                self.name + '.protein_domains.txt'
            ),
            'w'
        )
        fout.write(
            ','.join(['{}']*11).format(
                "5\'_gene",
                "3\'_gene",
                "5\'_transcript",
                "3\'_transcript",
                "5\'_strand",
                "3\'_strand",
                "Domain_ID",
                "Domain_name",
                "Domain_description",
                "Protein_start",
                "Protein_end"
            ) + '\n')

        for name, transcript in list(self.transcripts.items()):
            for domain in transcript.domains['fusion']:
                fout.write(
                    ','.join(['{}']*11).format(
                        transcript.gene5prime.gene.gene_name,
                        transcript.gene3prime.gene.gene_name,
                        transcript.transcript1.id,
                        transcript.transcript2.id,
                        transcript.transcript1.strand,
                        transcript.transcript2.strand,
                        domain[0],
                        domain[1],
                        domain[2],
                        domain[3],
                        domain[4]
                    ) + '\n')
        fout.close()

        # write tables of exon structure

        fout = open(
            os.path.join(
                out_dir,
                self.name + '.exons.txt'
            ),
            'w'
        )
        fout.write(
            ','.join(['{}']*11).format(
                "5\'_gene",
                "3\'_gene",
                "5\'_transcript",
                "3\'_transcript",
                "5\'_strand",
                "3\'_strand",
                "exon_gene_source",
                "exon_number",
                "exon_chr",
                "exon_start",
                "exon_end"
            ) + '\n')

        for name, transcript in list(self.transcripts.items()):
            for exon in transcript.gene5prime_exon_intervals:
                fout.write(
                    ','.join(['{}']*11).format(
                        transcript.gene5prime.gene.gene_name,
                        transcript.gene3prime.gene.gene_name,
                        transcript.transcript1.id,
                        transcript.transcript2.id,
                        transcript.transcript1.strand,
                        transcript.transcript2.strand,
                        '\'5 gene',
                        exon[2],
                        transcript.transcript1.contig,
                        exon[0],
                        exon[1]
                    ) + '\n')

            for exon in transcript.gene3prime_exon_intervals:
                fout.write(
                    ','.join(['{}']*11).format(
                        transcript.gene5prime.gene.gene_name,
                        transcript.gene3prime.gene.gene_name,
                        transcript.transcript1.id,
                        transcript.transcript2.id,
                        transcript.transcript1.strand,
                        transcript.transcript2.strand,
                        '\'3 gene',
                        exon[2],
                        transcript.transcript2.contig,
                        exon[0],
                        exon[1]
                    ) + '\n')

        fout.close()


class FusionTranscript(object):
    """
    Generates the information needed for the gene fusion transctips
    """

    def __init__(
            self,
            transcript1=None, transcript2=None,
            gene5prime=None, gene3prime=None,
            db=None, protein_databases=None,
            predict_effect=True):

        self.transcript1 = transcript1
        self.transcript2 = transcript2
        self.gene5prime = gene5prime
        self.gene3prime = gene3prime
        self.protein_databases = protein_databases
        self.db = db

        self.name = self.transcript1.id + '-' + self.transcript2.id
        self.gene_names = self.gene5prime.gene.name + '-' + self.gene3prime.gene.name

        self.cds = None
        self.transcript_cds_junction_5prime = None
        self.transcript_cds_junction_3prime = None
        self.gene5prime_cds_intervals = []
        self.gene3prime_cds_intervals = []
        self.cds_3prime = None
        self.cds_3prime = None

        self.cdna = None
        self.cdna_5prime = None
        self.cdna_3prime = None
        self.transcript_cdna_junction_5prime = None
        self.transcript_cdna_junction_3prime = None
        self.gene5prime_exon_intervals = []
        self.gene3prime_exon_intervals = []

        self.protein = None
        self.protein_length = None
        self.molecular_weight = None
        self.domains = {
            transcript1.id: [],
            transcript2.id: [],
            'fusion': []
        }

        self.transcript_protein_junction_5prime = None
        self.transcript_protein_junction_3prime = None

        self.effect_5prime = 'exon'
        self.effect_3prime = 'exon'
        self.effect = ''
        self.has_coding_potential = False

        if predict_effect:
            self.predict_effect()

    def _annotate(self):
        """
        Annotate the gene fusion's protein using the protein annotaiton
        from its two genes
        """

        fusion_domains = []
        gene5prime_domains = []
        gene3prime_domains = []

        tmp_domains = []

        # fetch the translation ids

        sqlite3_command = "SELECT * FROM " + self.db.build + "_transcript WHERE transcript_stable_id==\"" + self.transcript1.id + "\""
        self.db.logger.debug('SQLite - ' + sqlite3_command)
        self.db.sqlite3_cursor.execute(
            sqlite3_command
        )
        gene5prime_translation_id = self.db.sqlite3_cursor.fetchall()[0][3]

        sqlite3_command = "SELECT * FROM " + self.db.build + "_transcript WHERE transcript_stable_id==\"" + self.transcript2.id + "\""
        self.db.logger.debug('SQLite - ' + sqlite3_command)
        self.db.sqlite3_cursor.execute(
            sqlite3_command
        )
        gene3prime_translation_id = self.db.sqlite3_cursor.fetchall()[0][3]

        for protein_database in self.protein_databases:

            # fetch protein annotation

            sqlite3_command = "SELECT * FROM " + self.db.build + "_" + protein_database + " WHERE translation_id==\"" + gene5prime_translation_id + "\""
            self.db.logger.debug('SQLite - ' + sqlite3_command)
            self.db.sqlite3_cursor.execute(
                sqlite3_command
            )
            tmp_domains += [list(x) for x in self.db.sqlite3_cursor.fetchall()]

        for d in tmp_domains:

            pfeature_ID = d[2]
            pfeature_name = d[6]
            pfeature_description = d[5]
            pfeature_start = int(d[3])
            pfeature_end = int(d[4])

            gene5prime_domains.append([
                pfeature_ID,
                pfeature_name,
                pfeature_description,
                pfeature_start,
                pfeature_end
            ])

            try:
                if self.transcript_protein_junction_5prime < pfeature_start:
                    continue
                elif self.transcript_protein_junction_5prime >= pfeature_end:

                    fusion_domains.append([
                        pfeature_ID,
                        pfeature_name,
                        pfeature_description,
                        pfeature_start,
                        pfeature_end
                    ])

                elif (self.transcript_protein_junction_5prime - pfeature_start) >= MIN_DOMAIN_LENGTH:

                    pfeature_end = self.transcript_protein_junction_5prime

                    fusion_domains.append([
                        pfeature_ID,
                        pfeature_name,
                        pfeature_description,
                        pfeature_start,
                        pfeature_end
                    ])
            except:
                import pdb; pdb.set_trace()

        # only find the 3' gene partner's domains if the fusion is in-frame

        if self.effect != 'out-of-frame':

            tmp_domains = []

            for database in self.protein_databases:

                sqlite3_command = "SELECT * FROM " + self.db.build + "_" + database + " WHERE translation_id==\"" + gene3prime_translation_id + "\""
                self.db.logger.debug('SQLite - ' + sqlite3_command)
                self.db.sqlite3_cursor.execute(
                    sqlite3_command
                )
                tmp_domains += [list(x) for x in self.db.sqlite3_cursor.fetchall()]

            for d in tmp_domains:

                pfeature_ID = d[2]
                pfeature_name = d[6]
                pfeature_description = d[5]
                pfeature_start = int(d[3])
                pfeature_end = int(d[4])

                gene3prime_domains.append([
                    pfeature_ID,
                    pfeature_name,
                    pfeature_description,
                    pfeature_start,
                    pfeature_end
                ])

                if self.transcript_protein_junction_3prime > pfeature_end:
                    continue
                elif self.transcript_protein_junction_3prime <= pfeature_start:

                    pfeature_start = (pfeature_start-self.transcript_protein_junction_3prime) + self.transcript_protein_junction_5prime
                    pfeature_end = (pfeature_end-self.transcript_protein_junction_3prime) + self.transcript_protein_junction_5prime

                    fusion_domains.append([
                        pfeature_ID,
                        pfeature_name,
                        pfeature_description,
                        pfeature_start,
                        pfeature_end
                    ])

                else:

                    pfeature_start = self.transcript_protein_junction_5prime
                    pfeature_end = (pfeature_end-self.transcript_protein_junction_3prime) + self.transcript_protein_junction_5prime

                    fusion_domains.append([
                        pfeature_ID,
                        pfeature_name,
                        pfeature_description,
                        pfeature_start,
                        pfeature_end
                    ])

                if pfeature_end < pfeature_start:
                    import pdb; pdb.set_trace()

        self.domains['fusion'] = fusion_domains
        self.domains[self.transcript1.id] = gene5prime_domains
        self.domains[self.transcript2.id] = gene3prime_domains

    def _fetch_protein(self):
        """
        Predict the potential protain amino acid sequence
        """

        self.transcript_protein_junction_5prime = 0
        self.transcript_protein_junction_3prime = 0

        self.protein_names = self.transcript1.protein_id + '-' + self.transcript2.protein_id

        self.transcript_protein_junction_5prime = int(self.transcript_cds_junction_5prime/3.)

        if (len(self.cds_5prime)/3.).is_integer() and (len(self.cds_3prime)/3.).is_integer():
            self.effect='in-frame'
            self.transcript_protein_junction_3prime = int(self.transcript_cds_junction_3prime/3.)
        elif round((len(self.cds_5prime)/3. % 1) + (len(self.cds_3prime)/3. % 1),2) == 1.0:
            self.effect = 'in-frame (with mutation)'
            self.transcript_protein_junction_3prime = int(self.transcript_cds_junction_3prime/3.)
        else:
            self.effect = 'out-of-frame'

        # check if CDS's length is multiple of 3, if not then print warning

        if (len(self.cds.seq) % 3) !=0:
            self.db.logger.warn(
                'Length of fusion isoform CDS is not a multiple of 3!')

        # translate CDS into protein and remove everything after the stop codon

        if self.effect == 'out-of-frame':

            # trim the CDS sequence if fusion is out-of-frame

            seq = self.cds.seq[0:3*int(len(self.cds.seq)/3)]

            protein_seq = seq.translate()
            protein_seq = protein_seq[0:protein_seq.find('*')]
        else:
            protein_seq = self.cds.seq.translate()
            protein_seq = protein_seq[0:protein_seq.find('*')]

        # predict molecular weight

        self.molecular_weight = SeqUtils.molecular_weight(protein_seq)/1000.

        # convert to a sequence record

        self.protein_length = len(protein_seq)

        self.protein = SeqRecord.SeqRecord(
            protein_seq,
            id=self.protein_names,
            name=self.protein_names,
            description=("length: {}, kD: {}, transcripts: {}, strands: {}/{}, "
                         "genes: {}, effect: {}").format(
                            self.protein_length,
                            self.molecular_weight,
                            self.name,
                            self.transcript1.strand,
                            self.transcript2.strand,
                            self.gene_names,
                            self.effect
                         ))

    def _fetch_transcript_cds(self):
        """
        Predict the potential nucleotide sequence
        """

        self.transcript_cds_junction_5prime = 0
        self.transcript_cds_junction_3prime = 0

        #5prime transcript

        if self.transcript1.strand=="+":
            for cds in self.transcript1.coding_sequence_position_ranges:
                if self.gene5prime.junction >= cds[1]:
                    self.transcript_cds_junction_5prime += (cds[1] - cds[0] + 1)
                elif self.gene5prime.junction <= cds[0]:
                    break
                else:
                    self.transcript_cds_junction_5prime += (self.gene5prime.junction - cds[0] + 1)
                    break
        else:
            for cds in self.transcript1.coding_sequence_position_ranges:
                if self.gene5prime.junction <= cds[0]:
                    self.transcript_cds_junction_5prime += (cds[1] - cds[0] + 1)
                elif self.gene5prime.junction >= cds[1]:
                    break
                else:
                    self.transcript_cds_junction_5prime += (cds[1] - self.gene5prime.junction + 1)

        self.cds_5prime = self.transcript1.coding_sequence[0:self.transcript_cds_junction_5prime]

        if self.transcript2.strand=="+":
            for cds in self.transcript2.coding_sequence_position_ranges:
                if self.gene3prime.junction >= cds[1]:
                    self.transcript_cds_junction_3prime += (cds[1] - cds[0] + 1)
                elif self.gene3prime.junction <= cds[0]:
                    break
                else:
                    self.transcript_cds_junction_3prime += (self.gene3prime.junction - cds[0])
        else:
            for cds in self.transcript2.coding_sequence_position_ranges:
                if self.gene3prime.junction <= cds[0]:
                    self.transcript_cds_junction_3prime += (cds[1] - cds[0] + 1)
                elif self.gene3prime.junction >= cds[1]:
                    break
                else:
                    self.transcript_cds_junction_3prime += (cds[1] - self.gene3prime.junction)

        self.cds_3prime = self.transcript2.coding_sequence[self.transcript_cds_junction_3prime::]

        #create a sequence record

        if self.cds_5prime is not None and self.cds_3prime is not None:
            seq = self.cds_5prime + self.cds_3prime

        self.cds = SeqRecord.SeqRecord(
            Seq.Seq(seq,generic_dna),
            id=self.name,
            name=self.name,
            description="length: {}, genes: {}/{}, strands: {}/{}".format(
                len(self.cds_5prime + self.cds_3prime),
                self.transcript1.gene.name,
                self.transcript2.gene.name,
                self.transcript1.gene.strand,
                self.transcript2.gene.strand
            )
        )

    def _fetch_transcript_cdna(self):
        """
        Predict the potential nucleotide sequence
        """

        self.transcript_cdna_junction_5prime = 0
        self.transcript_cdna_junction_3prime = 0

        # get the 5prime transcript sequence and determine if junction is
        # within intron

        n = 0

        n_max = len(self.transcript1.exon_intervals)
        exons = self.transcript1.exon_intervals

        if self.transcript1.strand == "+":

            exon_count = 0

            for exon in exons:

                exon_count += 1

                # is in intron?

                if n == 0 and self.gene5prime.junction < exon[0]:
                    self.effect_5prime = 'intron (before cds)'

                if n == n_max and self.gene5prime.junction > exon[1]:
                    self.effect_5prime = 'intron (after cds)'
                elif self.gene5prime.junction > exon[1] and self.gene5prime.junction < exons[n+1][0]:
                    self.effect_5prime = 'intron'

                n += 1

                # get sequence

                if self.gene5prime.junction >= exon[1]:
                    self.transcript_cdna_junction_5prime += (exon[1] - exon[0] + 1)
                    self.gene5prime_exon_intervals.append([
                        exon[0],
                        exon[1],
                        exon_count
                    ])
                elif self.gene5prime.junction <= exon[0]:
                    break
                else:
                    self.transcript_cdna_junction_5prime += (self.gene5prime.junction - exon[0] + 1)
                    self.gene5prime_exon_intervals.append([
                        exon[0],
                        self.gene5prime.junction,
                        exon_count
                    ])
                    break
        else:

            exon_count = 0

            for exon in exons:

                exon_count += 1

                # is in intron?

                if n == 0 and self.gene5prime.junction > exon[1]:
                    self.effect_5prime = 'intron (before cds)'

                if n == n_max and self.gene5prime.junction < exon[0]:
                    self.effect_5prime = 'intron (after cds)'
                elif self.gene5prime.junction < exon[0] and self.gene5prime.junction > exons[n+1][1]:
                    self.effect_5prime = 'intron'

                n += 1

                # get sequence

                if self.gene5prime.junction <= exon[0]:
                    self.transcript_cdna_junction_5prime += (exon[1] - exon[0] + 1)
                    self.gene5prime_exon_intervals.append([
                        exon[0],
                        exon[1],
                        exon_count
                    ])
                elif self.gene5prime.junction >= exon[1]:
                    break
                else:
                    self.transcript_cdna_junction_5prime += (exon[1] - self.gene5prime.junction + 1)
                    self.gene5prime_exon_intervals.append([
                        self.gene5prime.junction,
                        exon[1],
                        exon_count
                    ])

        try:
            self.cdna_5prime = self.transcript1.sequence[0:self.transcript_cdna_junction_5prime]
        except TypeError:
            self.db.logger.warn('No cDNA sequence available for %s! ' \
                'Will not print cDNA sequence for the %s fusion. ' \
                'You might be working with an outdated pyensembl. ' \
                'Update the package and rerun \'pyensembl install\'' %
                (
                    str(self.gene5prime.gene.name),
                    self.gene_names
                )
            )

        # get the 3prime transcript sequence and determine if junction is
        # within intron

        n = 0
        n_max = len(self.transcript2.exon_intervals)
        exons = self.transcript2.exon_intervals

        if self.transcript2.strand == "+":

            exon_count = 0

            # get the exons in the fusion

            for exon in exons:
                exon_count += 1
                if self.gene3prime.junction >= exon[1]:
                    continue
                elif self.gene3prime.junction <= exon[0]:
                    self.gene3prime_exon_intervals.append([
                        exon[0],
                        exon[1],
                        exon_count
                    ])
                else:
                    self.gene3prime_exon_intervals.append([
                        self.gene3prime.junction,
                        exon[1],
                        exon_count
                    ])

            for exon in exons:

                # is in intron?

                if n == 0 and self.gene3prime.junction < exon[0]:
                    self.effect_3prime = 'intron (before cds)'

                if n == n_max and self.gene3prime.junction > exon[1]:
                    self.effect_3prime = 'intron (after cds)'
                elif self.gene3prime.junction > exon[1] and self.gene3prime.junction < exons[n+1][0]:
                    self.effect_3prime = 'intron'

                n += 1

                # get sequence

                if self.gene3prime.junction >= exon[1]:
                    self.transcript_cdna_junction_3prime += (exon[1] - exon[0] + 1)
                elif self.gene3prime.junction <= exon[0]:
                    break
                else:
                    self.transcript_cdna_junction_3prime += (self.gene3prime.junction - exon[0])
        else:

            exon_count = 0

            # get the exons in the fusion

            for exon in exons:
                exon_count += 1
                if self.gene3prime.junction <= exon[0]:
                    continue
                elif self.gene3prime.junction >= exon[1]:
                    self.gene3prime_exon_intervals.append([
                        exon[0],
                        exon[1],
                        exon_count
                    ])
                else:
                    self.gene3prime_exon_intervals.append([
                        exon[0],
                        self.gene3prime.junction,
                        exon_count
                    ])

            for exon in exons:

                # is in intron?

                if n == 0 and self.gene3prime.junction > exon[1]:
                    self.effect_3prime = 'intron (before cds)'

                if n == n_max and self.gene3prime.junction < exon[0]:
                    self.effect_3prime = 'intron (after cds)'
                elif self.gene3prime.junction < exon[0] and self.gene3prime.junction > exons[n+1][1]:
                    self.effect_3prime = 'intron'

                n += 1

                # get sequence

                if self.gene3prime.junction <= exon[0]:
                    self.transcript_cdna_junction_3prime += (exon[1] - exon[0] + 1)
                elif self.gene3prime.junction >= exon[1]:
                    break
                else:
                    self.transcript_cdna_junction_3prime += (exon[1] - self.gene3prime.junction)

        if self.cdna_5prime is not None:
            try:
                self.cdna_3prime = self.transcript2.sequence[
                    self.transcript_cdna_junction_3prime::
                ]
            except TypeError:
                self.db.logger.warn('No cDNA sequence available for %s! ' \
                    'Will not print cDNA sequence for the %s fusion. ' \
                    'You might be working with an outdated pyensembl. ' \
                    'Update the package and rerun \'pyensembl install\'' %
                    (
                        str(self.gene3prime.gene.name),
                        self.gene_names
                    )
                )

        if self.cdna_5prime is not None and self.cdna_3prime is not None:
            seq = self.cdna_5prime + self.cdna_3prime
            seq_length = str(len(seq))
        else:
            seq = ''
            seq_length = 'NA'

        self.cdna = SeqRecord.SeqRecord(
            Seq.Seq(seq, generic_dna),
            id=self.name,
            name=self.name,
            description="length=" + seq_length
        )

        # find out if the junction on either gene is with in the 5' or 3' UTR,
        # or if it exactly at the beginning or end of the UTR

        # the 5' gene

        if self.effect_5prime.find('intron') == -1:

            if self.transcript1.complete and \
                    self.transcript_cdna_junction_5prime < len(self.transcript1.five_prime_utr_sequence):
                self.effect_5prime='5UTR'
            elif self.transcript1.complete and \
                    self.transcript_cdna_junction_5prime == len(self.transcript1.five_prime_utr_sequence):
                self.effect_5prime='5UTR (end)'

            if self.transcript1.complete and \
                    (len(self.transcript1)-self.transcript_cdna_junction_5prime) < len(self.transcript1.three_prime_utr_sequence):
                self.effect_5prime='3UTR'
            elif self.transcript1.complete and \
                    (len(self.transcript1)-self.transcript_cdna_junction_5prime) == len(self.transcript1.three_prime_utr_sequence):
                self.effect_5prime='3UTR (start)'

        #the 3' gene

        if self.effect_3prime.find('intron')==-1:

            if self.transcript2.complete and \
                    self.transcript_cdna_junction_3prime < len(self.transcript2.five_prime_utr_sequence):
                self.effect_3prime='5UTR'
            elif self.transcript2.complete and \
                    self.transcript_cdna_junction_3prime == len(self.transcript2.five_prime_utr_sequence):
                self.effect_3prime='5UTR (end)'

            if self.transcript2.complete and \
                    (len(self.transcript2)-self.transcript_cdna_junction_3prime) < len(self.transcript2.three_prime_utr_sequence):
                self.effect_3prime='3UTR'
            elif self.transcript2.complete and \
                    (len(self.transcript2)-self.transcript_cdna_junction_3prime) == len(self.transcript2.three_prime_utr_sequence):
                self.effect_3prime='3UTR (start)'

        self.effect = self.effect_5prime + '-' + self.effect_3prime

    # def _check_if_in_intron(self):

    def predict_effect(self):
        """
        For all gene isoform combinations predict the effect of the fusion
        (e.g. in-frame, out-of-frame, etc...). Then if it has coding potential
        then annotate with domains
        """

        # check if within exon or UTR

        self._fetch_transcript_cdna()

        # check if within CDS and if it occurs at the very beginning or end of CDS

        if self.transcript1.complete and (self.effect_5prime.find('UTR')==-1):

            if self.effect_5prime.find('intron') != -1:

                # if in intron, is it between CDS regions?

                n = 0
                n_max = len(self.transcript1.coding_sequence_position_ranges)-1

                for cds in self.transcript1.coding_sequence_position_ranges:

                    if self.transcript1.strand == "+":
                        if n == 0 and self.gene5prime.junction < cds[0]:
                            self.effect_5prime = 'intron (before cds)'
                            break
                        elif n == n_max and self.gene5prime.junction > cds[1]:
                            self.effect_5prime = 'intron (after cds)'
                            break
                    else:
                        if n == 0 and self.gene5prime.junction > cds[1]:
                            self.effect_5prime = 'intron (before cds)'
                            break
                        elif n == n_max and self.gene5prime.junction < cds[0]:
                            self.effect_5prime = 'intron (after cds)'
                            break


                    if (n!=0) and \
                            (n!=n_max) and \
                            (
                                self.gene5prime.junction < cds[0] or \
                                self.gene5prime.junction > self.transcript1.coding_sequence_position_ranges[n+1][1]
                            ):

                        self.effect_5prime='intron (cds)'
                        break

                    n += 1
            else:
                for cds in self.transcript1.coding_sequence_position_ranges:

                    if self.gene5prime.junction >= cds[0] and self.gene5prime.junction <= cds[1]:
                        self.effect_5prime = 'CDS'
                        break

                if self.transcript1.strand == "+":
                    if self.gene5prime.junction==self.transcript1.coding_sequence_position_ranges[0][0]:
                        self.effect_5prime = 'CDS (start)'
                    elif self.gene5prime.junction==self.transcript1.coding_sequence_position_ranges[-1][1]:
                        self.effect_5prime = 'CDS (end)'
                else:
                    if self.gene5prime.junction==self.transcript1.coding_sequence_position_ranges[0][1]:
                        self.effect_5prime = 'CDS (start)'
                    elif self.gene5prime.junction==self.transcript1.coding_sequence_position_ranges[-1][0]:
                        self.effect_5prime = 'CDS (end)'

        if self.transcript2.complete and (self.effect_3prime.find('UTR') == -1):

            if self.effect_3prime.find('intron') != -1:

                #if in intron, is it between CDS regions?

                n = 0
                n_max = len(self.transcript2.coding_sequence_position_ranges)-1

                for cds in self.transcript2.coding_sequence_position_ranges:

                    if self.transcript2.strand=="+":
                        if n==0 and self.gene3prime.junction < cds[0]:
                            self.effect_3prime='intron (before cds)'
                            break
                        elif n==n_max and self.gene3prime.junction > cds[1]:
                            self.effect_3prime='intron (after cds)'
                            break
                    else:
                        if n == 0 and self.gene3prime.junction > cds[1]:
                            self.effect_3prime = 'intron (before cds)'
                            break
                        elif n == n_max and self.gene3prime.junction < cds[0]:
                            self.effect_3prime = 'intron (after cds)'
                            break

                    if (n != 0) and \
                            (n != n_max) and \
                            (
                                self.gene3prime.junction < cds[0] or \
                                self.gene3prime.junction > self.transcript2.coding_sequence_position_ranges[n+1][1]
                            ):

                        self.effect_3prime = 'intron (cds)'
                        break

                    n += 1

            else:
                for cds in self.transcript2.coding_sequence_position_ranges:
                    if self.gene3prime.junction >= cds[0] and self.gene3prime.junction <= cds[1]:
                        self.effect_3prime='CDS'
                        break

                if self.transcript2.strand == "+":
                    if self.gene3prime.junction==self.transcript2.coding_sequence_position_ranges[0][0]:
                        self.effect_3prime = 'CDS (start)'
                    elif self.gene3prime.junction==self.transcript2.coding_sequence_position_ranges[-1][1]:
                        self.effect_3prime = 'CDS (end)'
                else:
                    if self.gene3prime.junction==self.transcript2.coding_sequence_position_ranges[0][1]:
                        self.effect_3prime = 'CDS (start)'
                    elif self.gene3prime.junction==self.transcript2.coding_sequence_position_ranges[-1][0]:
                        self.effect_3prime = 'CDS (end)'

        # get if has coding potential

        self.has_coding_potential = utils.CODING_COMBINATIONS[(self.effect_5prime,self.effect_3prime)]['protein_coding_potential']

        # check if they don't have stop and/or start codons

        reasons = []

        if not self.transcript1.contains_start_codon:
            self.has_coding_potential=False
            reasons.append("no known 5' transcript start codon")
        if not self.transcript1.contains_stop_codon:
            self.has_coding_potential=False
            reasons.append("no known 5' transcript stop codon")

        if not self.transcript2.contains_start_codon:
            self.has_coding_potential=False
            reasons.append("no known 3' transcript start codon")
        if not self.transcript2.contains_stop_codon:
            self.has_coding_potential=False
            reasons.append("no known 3' transcript stop codon")

        # append information to cdna fasta headers

        self.cdna.description += "; locations: {}/{};".format(
            self.effect_5prime, self.effect_3prime)
        self.cdna.description += " strands: {}/{};".format(
            self.transcript1.strand, self.transcript2.strand)
        self.cdna.description += " Has protein coding potential: {};".format(
            self.has_coding_potential)

        if not self.has_coding_potential:
            self.cdna.description += " Reason: {}".format(', '.join(reasons))

        # if the fusion transcript has coding potential then
        # fetch its CDS and protein sequences and annotate it

        if self.has_coding_potential:
            self._fetch_transcript_cds()
            self._fetch_protein()
            self._annotate()
