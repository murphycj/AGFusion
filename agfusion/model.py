import itertools
import os

from agfusion import utils, exceptions
import pandas
from PIL import *
import pyensembl
from Bio import Seq, SeqIO, SeqRecord, SeqUtils
from Bio.Alphabet import generic_dna,generic_protein

class Model(object):
    def __new__(cls, *args, **kwargs):
        if cls is Model:
            raise TypeError("Model shoud not be instantiated")
        return object.__new__(cls, *args, **kwargs)

        self.is_fusion = False

        self.domains={
            i:{''} for i,j,k in utils.PROTEIN_DOMAIN
        }
        self.transmembrane={
            'transmembrane_domain':{
                'transmembrane_domain_start':0,
                'transmembrane_domain_end':0
            }
        }
        self.complexity={
            'low_complexity':{
                'low_complexity_start':0,
                'low_complexity_end':0
            }
        }

    def _does_dir_exist(self,out_dir):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    def save_image(self,out_dir,length_normalize=None,dpi=90,file_type='png',fontsize=12):

        import matplotlib.pyplot as plt
        import matplotlib.patches as patches

        self._does_dir_exist(out_dir)

        assert file_type in ['png','pdf','jpeg'], 'provided wrong file type'

        for name, transcript in self.transcripts.items():

            if not transcript.has_coding_potential:
                continue

            fig = plt.figure(figsize=(20,5))
            ax = fig.add_subplot(111)

            #main protein body

            ax.add_patch(
                patches.Rectangle(
                    (0.05, 0.45),
                    0.9,
                    0.1,
                    fill=False
                )
            )
            ax.text(
                0.5,
                0.9,
                name,
                horizontalalignment='center',
                fontsize=fontsize
            )

            #plot domains

            for domain in transcript.domains['pfam']:

                if length_normalize is None:
                    length_normalize = transcript.protein_length
                else:
                    assert length_normalize >= transcript.protein_length, "length normalization should be ≥ protein length"

                offset=0.05

                domain_name = domain[1]
                domain_start = (int(domain[2])/float(length_normalize))*0.9 + offset
                domain_end = (int(domain[3])/float(length_normalize))*0.9 + offset
                domain_center = (domain_end-domain_start)/2. + domain_start

                ax.add_patch(
                    patches.Rectangle(
                        (
                            domain_start,
                            0.45
                        ),
                        domain_end-domain_start,
                        0.1
                    )
                )
                ax.text(
                    domain_center,
                    0.4,
                    domain_name,
                    horizontalalignment='center',
                    fontsize=fontsize
                )

            #add the junction

            try:
                ax.add_line(plt.Line2D(
                    (
                        (transcript.transcript_protein_junction_5prime/float(length_normalize))*0.9 + offset,
                        (transcript.transcript_protein_junction_5prime/float(length_normalize))*0.9 + offset
                    ),
                    (0.4,0.6),
                    color='black'
                    )
                )
            except:
                import pdb; pdb.set_trace()

            ax.axis('off')
            ax.set_xlim(0,1)
            ax.set_ylim(0,1)
            filename = out_dir + '/' + name + '.' + file_type
            fig.savefig(
                filename,
                dpi=dpi,
                bbox_inches='tight'
            )
            plt.close(fig)
            plt.clf()

            length_normalize = None


    def save_transcript_cdna(self,out_dir):

        self._does_dir_exist(out_dir)

        fout = open(out_dir + '/cdna.fa','w')

        for name, transcript in self.transcripts.items():

            if transcript.cds is not None:

                SeqIO.write(transcript.cdna,fout,"fasta")

        fout.close()

    def save_transcript_cds(self,out_dir):

        self._does_dir_exist(out_dir)

        fout = open(out_dir + '/cds.fa','w')

        for name, transcript in self.transcripts.items():

            if transcript.cds is not None:

                SeqIO.write(transcript.cds,fout,"fasta")

        fout.close()

    def save_proteins(self,out_dir):

        self._does_dir_exist(out_dir)

        fout = open(out_dir + '/protein.fa','w')

        for name, transcript in self.transcripts.items():

            if transcript.cds is not None:

                SeqIO.write(transcript.protein,fout,"fasta")

        fout.close()

    def save_annotations(self,out_file,annotation='pfam'):

        fout = open(out_file,'w')
        fout.write("Gene,transcript,domain,protein_start,protein_end\n")

        for name, transcript in self.transcripts.items():
            for domain in transcript.domains[annotation]:
                try:
                    fout.write(
                        self.name + ',' + \
                        name + ',' + \
                        domain[1] + ',' + \
                        str(domain[2]) + ',' + \
                        str(domain[3]) + '\n'
                    )
                except:
                    import pdb; pdb.set_trace()
        fout.close()

class Gene(Model):
    """
    Stores the necessary information to specify the architecture of either
    wild-type genes or fusion gene.
    """

    def __init__(self,gene=None,junction=0,db=None):

        self.gene = gene
        self.db=db

        self.junction=junction
        if self.junction < self.gene.start or self.junction > self.gene.end:
            raise exceptions.JunctionException()

class Fusion(Model):
    """
    Generates the information needed for the gene fusion
    """

    def __init__(self,gene5prime=None,gene3prime=None,db=None,transcripts_5prime=None,transcripts_3prime=None):
        self.gene5prime=gene5prime
        self.gene3prime=gene3prime
        self.name = self.gene5prime.gene.name + '-' + self.gene3prime.gene.name
        self.db=db
        self.is_fusion = True

        self.transcripts = {}

        for combo in list(itertools.product(
                        self.gene5prime.gene.transcripts,
                        self.gene3prime.gene.transcripts)):
            transcript1 = combo[0]
            transcript2 = combo[1]

            if transcripts_5prime is not None and transcript1.id not in transcripts_5prime:
                continue
            if transcripts_3prime is not None and transcript2.id not in transcripts_3prime:
                continue

            #skip if the junction is outside the range of either transcript

            if not transcript1.contains(transcript1.contig,self.gene5prime.junction,self.gene5prime.junction):
                continue

            if not transcript2.contains(transcript2.contig,self.gene3prime.junction,self.gene3prime.junction):
                continue

            name = transcript1.id+'-'+transcript2.id
            self.transcripts[name] = FusionTranscript(
                transcript1,
                transcript2,
                gene5prime,
                gene3prime,
                db=db
            )

class FusionTranscript():
    """
    Generates the information needed for the gene fusion transctips
    """

    def __init__(self,transcript1=None,transcript2=None,gene5prime=None,gene3prime=None,db=None):
        self.transcript1=transcript1
        self.transcript2=transcript2
        self.gene5prime=gene5prime
        self.gene3prime=gene3prime
        self.db=db

        self.name=self.transcript1.id + '-' + self.transcript2.id
        self.gene_names = self.gene5prime.gene.name + '-' + self.gene3prime.gene.name

        self.transcript_protein_molecular_weight = 0
        self.cds=None
        self.cdna=None
        self.protein=None
        self.protein_length=None
        self.molecular_weight=None
        self.domains={
            i:[] for i,j,k in utils.PROTEIN_DOMAIN
        }
        self.transmembrane={
            'transmembrane_domain':{
                'transmembrane_domain_start':0,
                'transmembrane_domain_end':0
            }
        }
        self.complexity={
            'low_complexity':{
                'low_complexity_start':0,
                'low_complexity_end':0
            }
        }

        self.transcript_cdna_junction_5prime = 0
        self.transcript_cdna_junction_3prime = 0
        self.transcript_cds_junction_5prime = 0
        self.transcript_cds_junction_3prime = 0
        self.transcript_protein_junction_5prime = 0
        self.transcript_protein_junction_3prime = 0

        self.effect_5prime='exon'
        self.effect_3prime='exon'
        self.effect=''

        self.has_start_codon_5prime=True
        self.has_start_codon_3prime=True
        self.has_stop_codon_5prime=True
        self.has_stop_codon_3prime=True

        self.has_coding_potential=True
        self.predict_effect()

    def _fetch_domain_name(self,name):
        self.db.c.execute(
            'SELECT * FROM PFAMMAP WHERE pfam_acc==\"' + \
            name + \
            '\"'
        )
        new_name = self.db.c.fetchall()
        if len(new_name)<1:
            #import pdb; pdb.set_trace()
            #assert len(new_name)==1, "more than one name returned for " + d[1]
            print 'No pfam name for %s' % name
            return name
        if len(new_name)>1:
            import pdb; pdb.set_trace()
        elif len(new_name)==1:
            return new_name[0][1]

    def _annotate(self):
        """
        Annotate the gene fusion's protein using the protein annotaiton
        from its two genes
        """

        self.db.c.execute(
            'SELECT * FROM pfam WHERE ensembl_transcript_id==\"' + \
            self.transcript1.id + \
            '\"'
        )

        fusion_domains=[]

        domains = map(lambda x: list(x),self.db.c.fetchall())

        for d in domains:
            if str(d[1])=='':
                continue

            d[0]=self.name
            d[1] = self._fetch_domain_name(d[1])
            d[2] = int(d[2])
            d[3] = int(d[3])

            try:
                if self.transcript_protein_junction_5prime < d[2]:
                    continue
                elif self.transcript_protein_junction_5prime >= d[3]:
                    fusion_domains.append(d)
                else:
                    d[3] = self.transcript_protein_junction_5prime
                    fusion_domains.append(d)
            except:
                import pdb; pdb.set_trace()

        self.db.c.execute(
            'SELECT * FROM pfam WHERE ensembl_transcript_id==\"' + \
            self.transcript2.id + \
            '\"'
        )

        domains = map(lambda x: list(x),self.db.c.fetchall())
        for d in domains:
            if str(d[1])=='':
                continue

            d[0]=self.name
            d[1] = self._fetch_domain_name(d[1])
            d[2] = int(d[2])
            d[3] = int(d[3])

            if self.transcript_protein_junction_3prime > d[3]:
                continue
            elif self.transcript_protein_junction_5prime <= d[2]:
                d[2]=(d[2]-self.transcript_protein_junction_3prime) + self.transcript_protein_junction_5prime
                d[3]=(d[3]-self.transcript_protein_junction_3prime) + self.transcript_protein_junction_5prime
                fusion_domains.append(d)
            else:
                d[2]=(d[2]-self.transcript_protein_junction_3prime) + self.transcript_protein_junction_5prime
                d[3]=(d[3]-self.transcript_protein_junction_3prime) + self.transcript_protein_junction_5prime
                fusion_domains.append(d)

            if d[3]<d[2]:
                import pdb; pdb.set_trace()


        self.domains['pfam'] = fusion_domains

    def _fetch_protein(self):
        """
        Predict the potential protain amino acid sequence
        """

        self.transcript_protein_junction_5prime = 0
        self.transcript_protein_junction_3prime = 0

        self.protein_names = self.transcript1.protein_id + '-' + self.transcript2.protein_id

        #TODO the fusion cds could be a multiple of three but the 5' and 3'
        # cds are not.

        if (len(self.cds_5prime)/3.).is_integer() and (len(self.cds_3prime)/3.).is_integer():
            self.effect='in-frame'
        elif (len(self.cds)/3.).is_integer():
            self.effect='in-frame'
        else:
            self.effect='out-of-frame'

        self.transcript_protein_junction_5prime = int(self.transcript_cds_junction_5prime/3.)
        #self.transcript_protein_junction_3prime = len(self.transcript2.protein_sequence) - int(self.transcript_cds_junction_3prime/3.)
        self.transcript_protein_junction_3prime = int(self.transcript_cds_junction_3prime/3.)

        protein_seq = self.cds.seq.translate()
        protein_seq = protein_seq[0:protein_seq.find('*')]

        self.molecular_weight = SeqUtils.molecular_weight(protein_seq)/1000.
        self.protein_length = len(protein_seq)

        self.protein = SeqRecord.SeqRecord(
            protein_seq,
            id=self.protein_names,
            name=self.protein_names,
            description="length=" + str(self.protein_length) + \
                ", kD: " + str(self.molecular_weight) + \
                ", transcripts: " + str(self.name) + ', ' + \
                ", genes: " + str(self.gene_names)
        )

        self._annotate()

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
                    self.transcript_cds_junction_3prime += (cds[1] - self.gene3prime.junction + 1)

        self.cds_3prime = self.transcript2.coding_sequence[self.transcript_cds_junction_3prime::]

        self.cds = SeqRecord.SeqRecord(
            Seq.Seq(self.cds_5prime + self.cds_3prime,generic_dna),
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
        transcript_seq=''

        #5prime transcript

        if self.transcript1.strand=="+":
            for exon in self.transcript1.exons:
                if self.gene5prime.junction >= exon.end:
                    self.transcript_cdna_junction_5prime += (exon.end - exon.start + 1)
                elif self.gene5prime.junction <= exon.start:
                    break
                else:
                    self.transcript_cdna_junction_5prime+=self.gene5prime.junction-exon.start
                    break
        else:
            for exon in self.transcript1.exons:
                if self.gene5prime.junction <= exon.start:
                    self.transcript_cdna_junction_5prime += (exon.end - exon.start + 1)
                elif self.gene5prime.junction >= exon.end:
                    break
                else:
                    self.transcript_cdna_junction_5prime += (exon.end - self.gene5prime.junction + 1)

        transcript_seq += self.transcript1.sequence[0:self.transcript_cdna_junction_5prime]

        if self.transcript2.strand=="+":
            for exon in self.transcript2.exons:
                if self.gene3prime.junction >= exon.end:
                    self.transcript_cdna_junction_3prime += (exon.end - exon.start + 1)
                elif self.gene3prime.junction <= exon.start:
                    break
                else:
                    self.transcript_cdna_junction_3prime += (self.gene3prime.junction - exon.start + 1)
        else:
            for exon in self.transcript2.exons:
                if self.gene3prime.junction <= exon.start:
                    self.transcript_cdna_junction_3prime += (exon.end - exon.start + 1)
                elif self.gene3prime.junction >= exon.end:
                    break
                else:
                    self.transcript_cdna_junction_3prime += (exon.end - self.gene3prime.junction + 1)

        transcript_seq += self.transcript2.sequence[self.transcript_cdna_junction_3prime::]

        self.cdna = SeqRecord.SeqRecord(
            Seq.Seq(transcript_seq,generic_dna),
            id=self.name,
            name=self.name,
            description="length=" + str(len(transcript_seq)) + \
                ", alternative transcript names: " + str(self.transcript1.name) + ', ' + \
                str(self.transcript2.name)
        )

        if self.transcript1.complete and self.transcript_cdna_junction_5prime < len(self.transcript1.five_prime_utr_sequence):
            self.effect_5prime='5UTR'

        if self.transcript2.complete and (len(self.transcript2)-self.transcript_cdna_junction_5prime) < len(self.transcript2.five_prime_utr_sequence):
            self.effect_5prime='3UTR'

    def predict_effect(self):
        """
        types of fusions: in-frame, out-of-frame, and any combination of
        UTR, intron, intergenic, exonic (no known CDS),
        CDS(not-reliable-start-or-end), CDS(truncated), CDS(complete)


        CDS
        exon 5' CDS
        exon 3' CDS
        5' utr
        3' utr
        intron
        intergenic

        """

        #check if within exon

        self._fetch_transcript_cdna()

        #check if within CDS

        try:
            for cds in self.transcript1.coding_sequence_position_ranges:
                if self.gene5prime.junction >= cds[0] and self.gene5prime.junction <= cds[1]:
                    self.effect_5prime='CDS'
                    break
        except:
            pass

        try:
            for cds in self.transcript2.coding_sequence_position_ranges:
                if self.gene3prime.junction >= cds[0] and self.gene3prime.junction <= cds[1]:
                    self.effect_3prime='CDS'
                    break
        except:
            pass

        #if within CDS, then get CDS and predict protein

        if self.effect_5prime=='CDS' and self.effect_3prime=='CDS':
            self.has_coding_potential=True

        #check if the transcripts are complete, and if not check if
        #they don't have stop and/or start codons

        if not self.transcript1.contains_start_codon:
            self.has_start_codon_5prime=False
            self.has_coding_potential=False
        if not self.transcript1.contains_stop_codon:
            self.has_stop_codon_5prime=False
            self.has_coding_potential=False

        if not self.transcript2.contains_start_codon:
            self.has_start_codon_3prime=False
            self.has_coding_potential=False
        if not self.transcript2.contains_stop_codon:
            self.stop_codon_3prime=False
            self.has_coding_potential=False

        if self.has_coding_potential:
            self._fetch_transcript_cds()
            self._fetch_protein()
