from agfusion import utils, exceptions, plot
import pandas
from Bio import Seq, SeqIO, SeqRecord, SeqUtils
from Bio.Alphabet import generic_dna, generic_protein

import itertools
import os

# matplotlib.rcParams['interactive'] = False

MIN_DOMAIN_LENGTH = 5


class _Gene():
    """
    Stores the necessary information to specify the architecture of either
    wild-type genes or fusion gene.
    """

    def __init__(self, gene=None, junction=0, pyensembl_data=None,
                 genome='', gene5prime=False, db=None):

        gene = gene.upper()

        self.domains = []
        self.transcript_ids = {}

        # search for ensembl gene id first

        if gene in pyensembl_data.gene_ids():

            self.gene = pyensembl_data.gene_by_id(gene)

        else:

            # search by gene symbol next

            if pyensembl_data.species.latin_name == 'mus_musculus':
                gene = gene.capitalize()

            if gene in pyensembl_data.gene_names():
                temp = pyensembl_data.genes_by_name(gene)

                if len(temp) > 1:

                    # if too many ensembl gene IDs returned

                    ids = map(lambda x: x.id, temp)

                    raise exceptions.TooManyGenesException(gene, ids, genome)
                else:
                    self.gene = temp[0]
            else:
                if gene5prime:
                    raise exceptions.GeneIDException5prime(gene)
                else:
                    raise exceptions.GeneIDException3prime(gene)

        self.junction = junction

        if self.junction < self.gene.start or self.junction > self.gene.end:
            if gene5prime:
                raise exceptions.JunctionException5prime()
            else:
                raise exceptions.JunctionException3prime()

        # fetch the entrez gene id and canonical transcript id

        sqlite3_command = "SELECT * FROM " + db.build + " WHERE stable_id==\"" + self.gene.gene_id + "\""
        db.logger.debug('SQLite - ' + sqlite3_command)
        db.sqlite3_cursor.execute(
            sqlite3_command
        )
        tmp = db.sqlite3_cursor.fetchall()

        self.entrez_id = tmp[0][2]
        self.gene_id = tmp[0][0]

        sqlite3_command = "SELECT * FROM " + db.build + "_transcript WHERE transcript_id==\"" + tmp[0][4] + "\""
        db.logger.debug('SQLite - ' + sqlite3_command)
        db.sqlite3_cursor.execute(
            sqlite3_command
        )
        tmp = db.sqlite3_cursor.fetchall()

        self.canonical_transcript_id = tmp[0][2]

        # fetch transcript ids

        sqlite3_command = "SELECT * FROM " + db.build + "_transcript WHERE gene_id==\"" + self.gene_id + "\""
        db.logger.debug('SQLite - ' + sqlite3_command)
        db.sqlite3_cursor.execute(
            sqlite3_command
        )
        tmp = db.sqlite3_cursor.fetchall()

        for transcript in tmp:
            self.transcript_ids[transcript[2]] = transcript[0]


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

        # get the reference genom

        self.gene5prime = _Gene(
            gene=gene5prime,
            junction=gene5primejunction,
            pyensembl_data=pyensembl_data,
            gene5prime=True,
            db=db
        )

        self.gene3prime = _Gene(
            gene=gene3prime,
            junction=gene3primejunction,
            pyensembl_data=pyensembl_data,
            gene5prime=False,
            db=db
        )

        self.name = self.gene5prime.gene.name + '-' + self.gene3prime.gene.name

        # fetch which are canonical transcripts

        gene5prime_canonical = ''
        gene3prime_canonical = ''

        # construct all the fusion transcript combinations

        self.transcripts = {}

        for combo in list(itertools.product(
                        self.gene5prime.gene.transcripts,
                        self.gene3prime.gene.transcripts)):
            transcript1 = combo[0]
            transcript2 = combo[1]

            if transcripts_5prime is not None \
                    and transcript1.id not in transcripts_5prime:

                continue

            if transcripts_3prime is not None \
                    and transcript2.id not in transcripts_3prime:

                continue

            if (not noncanonical) and \
                    (transcript1.id != self.gene5prime.canonical_transcript_id):
                continue

            if (not noncanonical) and \
                    (transcript2.id != self.gene3prime.canonical_transcript_id):
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

        for name, transcript in self.transcripts.items():

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
            no_domain_labels=True, plot_WT=False):
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

        for name, transcript in self.transcripts.items():

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
                transcript=transcript
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
                    transcript=transcript
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
                    transcript=transcript
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

        for name, transcript in self.transcripts.items():

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
        for name, transcript in self.transcripts.items():

            if transcript.cds is not None:
                n+=1

        if n == 0:
            print 'Warning: the ' + self.name + ' fusion does not produce any protein coding transcripts. No cds.fa file will be written'
            return

        fout = open(
            os.path.join(
                out_dir,
                self.name + '_cds.fa'
            ),
            'w'
        )

        for name, transcript in self.transcripts.items():

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
        for name, transcript in self.transcripts.items():

            if transcript.cds is not None:
                n += 1

        if n == 0:
            print 'Warning: the ' + self.name + ' fusion does not produce any protein coding transcripts. No proteins.fa file will be written'
            return

        fout = open(
            os.path.join(
                out_dir,
                self.name + '_protein.fa'
            ),
            'w'
        )

        for name, transcript in self.transcripts.items():

            if transcript.cds is not None:

                if middlestar:
                    temp = str(transcript.protein.seq)
                    temp = temp[:transcript.transcript_protein_junction_5prime] + '*' + temp[transcript.transcript_protein_junction_5prime:]
                    transcript.protein.seq = Seq.Seq(temp,generic_protein)

                SeqIO.write(transcript.protein,fout,"fasta")

        fout.close()

#    def save_annotations(self,out_file='.',annotation='pfam'):
#
#        fout = open(out_file,'w')
#        fout.write("Gene,transcript,domain,protein_start,protein_end\n")
#
#        for name, transcript in self.transcripts.items():
#            for domain in transcript.domains[annotation]:
#                try:
#                    fout.write(
#                        self.name + ',' + \
#                        name + ',' + \
#                        domain[1] + ',' + \
#                        str(domain[2]) + ',' + \
#                        str(domain[3]) + '\n'
#                    )
#                except:
#                    import pdb; pdb.set_trace()
#        fout.close()


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

        self.cdna = None
        self.cdna_5prime = ''
        self.cdna_3prime = ''
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
        gene5prime_translation_id = self.db.sqlite3_cursor.fetchall()[0][4]

        sqlite3_command = "SELECT * FROM " + self.db.build + "_transcript WHERE transcript_stable_id==\"" + self.transcript2.id + "\""
        self.db.logger.debug('SQLite - ' + sqlite3_command)
        self.db.sqlite3_cursor.execute(
            sqlite3_command
        )
        gene3prime_translation_id = self.db.sqlite3_cursor.fetchall()[0][4]

        for protein_database in self.protein_databases:

            # fetch protein annotation

            sqlite3_command = "SELECT * FROM " + self.db.build + "_" + protein_database + " WHERE translation_id==\"" + gene5prime_translation_id + "\""
            self.db.logger.debug('SQLite - ' + sqlite3_command)
            self.db.sqlite3_cursor.execute(
                sqlite3_command
            )
            tmp_domains += map(lambda x: list(x), self.db.sqlite3_cursor.fetchall())

        for d in tmp_domains:

            pfeature_ID = d[2]
            pfeature_description = d[5]
            pfeature_start = int(d[3])
            pfeature_end = int(d[4])

            gene5prime_domains.append([
                pfeature_ID,
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
                        pfeature_description,
                        pfeature_start,
                        pfeature_end
                    ])

                elif (self.transcript_protein_junction_5prime - pfeature_start) >= MIN_DOMAIN_LENGTH:

                    pfeature_end = self.transcript_protein_junction_5prime

                    fusion_domains.append([
                        pfeature_ID,
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
                tmp_domains += map(lambda x: list(x), self.db.sqlite3_cursor.fetchall())

            for d in tmp_domains:

                pfeature_ID = d[2]
                pfeature_description = d[5]
                pfeature_start = int(d[3])
                pfeature_end = int(d[4])

                gene3prime_domains.append([
                    pfeature_ID,
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
                        pfeature_description,
                        pfeature_start,
                        pfeature_end
                    ])

                else:

                    pfeature_start = self.transcript_protein_junction_5prime
                    pfeature_end = (pfeature_end-self.transcript_protein_junction_3prime) + self.transcript_protein_junction_5prime

                    fusion_domains.append([
                        pfeature_ID,
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

        #translate CDS into protein and remove everything after the stop codon

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
            description="length=" + str(self.protein_length) + \
                ", kD: " + str(self.molecular_weight) + \
                ", transcripts: " + str(self.name) + \
                ", genes: " + str(self.gene_names) + \
                ", effect: " + self.effect
        )

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
                    self.transcript_cds_junction_3prime += (self.gene3prime.junction - cds[0] + 1)
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

        seq = self.cds_5prime + self.cds_3prime

        self.cds = SeqRecord.SeqRecord(
            Seq.Seq(seq,generic_dna),
            id=self.name,
            name=self.name,
            description="length=" + str(len(self.cds_5prime + self.cds_3prime)) + \
                ", genes: " + \
                str(self.transcript1.gene.name) + '-' + \
                str(self.transcript2.gene.name)
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
        n_max = len(self.transcript1.exons)

        if self.transcript1.strand == "+":
            for exon in self.transcript1.exons:

                # is in intron?

                if n == 0 and self.gene5prime.junction < exon.start:
                    self.effect_5prime = 'intron (before cds)'

                if n == n_max and self.gene5prime.junction > exon.end:
                    self.effect_5prime = 'intron (after cds)'
                elif self.gene5prime.junction > exon.end and self.gene5prime.junction < self.transcript1.exons[n+1].start:
                    self.effect_5prime = 'intron'

                n += 1

                # get sequence

                if self.gene5prime.junction >= exon.end:
                    self.transcript_cdna_junction_5prime += (exon.end - exon.start + 1)
                    self.gene5prime_exon_intervals.append([
                        exon.start,
                        exon.end
                    ])
                elif self.gene5prime.junction <= exon.start:
                    break
                else:
                    self.transcript_cdna_junction_5prime += (self.gene5prime.junction - exon.start + 1)
                    self.gene5prime_exon_intervals.append([
                        exon.start,
                        self.gene5prime.junction
                    ])
                    break
        else:
            for exon in self.transcript1.exons:

                # is in intron?

                if n == 0 and self.gene5prime.junction > exon.end:
                    self.effect_5prime = 'intron (before cds)'

                if n == n_max and self.gene5prime.junction < exon.start:
                    self.effect_5prime = 'intron (after cds)'
                elif self.gene5prime.junction < exon.start and self.gene5prime.junction > self.transcript1.exons[n+1].end:
                    self.effect_5prime = 'intron'

                n += 1

                # get sequence

                if self.gene5prime.junction <= exon.start:
                    self.transcript_cdna_junction_5prime += (exon.end - exon.start + 1)
                    self.gene5prime_exon_intervals.append([
                        exon.start,
                        exon.end
                    ])
                elif self.gene5prime.junction >= exon.end:
                    break
                else:
                    self.transcript_cdna_junction_5prime += (exon.end - self.gene5prime.junction + 1)
                    self.gene5prime_exon_intervals.append([
                        self.gene5prime.junction,
                        exon.end
                    ])

        self.cdna_5prime = self.transcript1.sequence[0:self.transcript_cdna_junction_5prime]

        # get the 3prime transcript sequence and determine if junction is
        # within intron

        n = 0
        n_max = len(self.transcript2.exons)

        if self.transcript2.strand == "+":

            # get the exons in the fusion

            for exon in self.transcript2.exons:
                if self.gene3prime.junction >= exon.end:
                    continue
                elif self.gene3prime.junction <= exon.start:
                    self.gene3prime_exon_intervals.append([
                        exon.start,
                        exon.end
                    ])
                else:
                    self.gene3prime_exon_intervals.append([
                        self.gene3prime.junction,
                        exon.end
                    ])

            for exon in self.transcript2.exons:

                # is in intron?

                if n == 0 and self.gene3prime.junction < exon.start:
                    self.effect_3prime = 'intron (before cds)'

                if n == n_max and self.gene3prime.junction > exon.end:
                    self.effect_3prime = 'intron (after cds)'
                elif self.gene3prime.junction > exon.end and self.gene3prime.junction < self.transcript2.exons[n+1].start:
                    self.effect_3prime = 'intron'

                n += 1

                # get sequence

                if self.gene3prime.junction >= exon.end:
                    self.transcript_cdna_junction_3prime += (exon.end - exon.start + 1)
                elif self.gene3prime.junction <= exon.start:
                    break
                else:
                    self.transcript_cdna_junction_3prime += (self.gene3prime.junction - exon.start)
        else:

            # get the exons in the fusion

            for exon in self.transcript2.exons:
                if self.gene3prime.junction <= exon.start:
                    continue
                elif self.gene3prime.junction >= exon.end:
                    self.gene3prime_exon_intervals.append([
                        exon.start,
                        exon.end
                    ])
                else:
                    self.gene3prime_exon_intervals.append([
                        exon.start,
                        self.gene3prime.junction
                    ])

            for exon in self.transcript2.exons:

                # is in intron?

                if n == 0 and self.gene3prime.junction > exon.end:
                    self.effect_3prime = 'intron (before cds)'

                if n == n_max and self.gene3prime.junction < exon.start:
                    self.effect_3prime = 'intron (after cds)'
                elif self.gene3prime.junction < exon.start and self.gene3prime.junction > self.transcript2.exons[n+1].end:
                    self.effect_3prime = 'intron'

                n += 1

                # get sequence

                if self.gene3prime.junction <= exon.start:
                    self.transcript_cdna_junction_3prime += (exon.end - exon.start + 1)
                elif self.gene3prime.junction >= exon.end:
                    break
                else:
                    self.transcript_cdna_junction_3prime += (exon.end - self.gene3prime.junction)

        self.cdna_3prime = self.transcript2.sequence[
            self.transcript_cdna_junction_3prime::
        ]

        self.cdna = SeqRecord.SeqRecord(
            Seq.Seq(self.cdna_5prime+self.cdna_3prime, generic_dna),
            id=self.name,
            name=self.name,
            description="length=" +
                        str(len(self.cdna_5prime+self.cdna_3prime)) +
                        ', ' +
                        str(self.transcript2.name)
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
                n_max = len(self.transcript1.exons)

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
                n_max = len(self.transcript2.exons)

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

        # check if not check if they don't have stop and/or start codons

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

        self.cdna.description += "; 5' location: " + self.effect_5prime + \
            "; 3' location: " + self.effect_3prime + \
            "; Has protein coding potential: " + str(self.has_coding_potential)

        if not self.has_coding_potential:
            self.cdna.description += "; Reason: " + ', '.join(reasons)

        # if the fusion transcript has coding potential then
        # fetch its CDS and protein sequences and annotate it

        if self.has_coding_potential:
            self._fetch_transcript_cds()
            self._fetch_protein()
            self._annotate()
