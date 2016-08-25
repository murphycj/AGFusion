import itertools

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

    def save_image(self,file_prefix):
        pass

class Gene(Model):
    """
    Stores the necessary information to specify the architecture of either
    wild-type genes or fusion gene.
    """

    def __init__(self,gene=None,junction=0,reference='',db=None):

        self.gene = gene
        self.reference=reference
        self.db=db

        self.junction=junction
        if self.junction < self.gene.start or self.junction > self.gene.end:
            raise exceptions.JunctionException()

class Fusion(Model):
    """
    Generates the information needed for the gene fusion
    """

    def __init__(self,gene5prime=None,gene3prime=None,db=None):
        self.gene5prime=gene5prime
        self.gene3prime=gene3prime
        self.db=db

        self.transcripts = {}

        for combo in list(itertools.product(
                        self.gene5prime.gene.transcripts,
                        self.gene3prime.gene.transcripts)):
            transcript1 = combo[0]
            transcript2 = combo[1]
            name=transcript1.id+'-'+transcript2.id
            self.transcripts[name] = FusionTranscript(
                transcript1,
                transcript2,
                gene5prime,
                gene3prime,
                db=db
            )

    def save_transcript_cdna(self,out_file):

        fout = open(out_file,'w')

        for name, fusion_transcript in self.transcripts.items():

            fusion_transcript = self.transcripts[name]

            fusion_transcript.fetch_transcript_cdna()

            SeqIO.write(fusion_transcript.cdna,fout,"fasta")

        fout.close()

    def save_transcript_cds(self,out_file):

        fout = open(out_file,'w')

        for name, fusion_transcript in self.transcripts.items():

            if not fusion_transcript.transcript1.complete or not fusion_transcript.transcript2.complete:
                continue

            fusion_transcript = self.transcripts[name]

            fusion_transcript.fetch_transcript_cds()

            SeqIO.write(fusion_transcript.cds,fout,"fasta")

        fout.close()

    def save_proteins(self,out_file):

        fout = open(out_file,'w')

        for name, fusion_transcript in self.transcripts.items():

            if not fusion_transcript.transcript1.complete or not fusion_transcript.transcript2.complete:
                continue

            fusion_transcript = self.transcripts[name]

            if fusion_transcript.cds is not None:
                fusion_transcript.fetch_transcript_cds()

            fusion_transcript.fetch_protein_seq()

            SeqIO.write(fusion_transcript.protein,fout,"fasta")

        fout.close()

class FusionTranscript(Model):
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

        self.transcript_cds = None
        self.transcript_cdna = None
        self.transcript_protein = None
        self.transcript_protein_molecular_weight = 0
        self.cds=None
        self.cdna=None
        self.protein=None
        self.domains={
            i:{} for i,j,k in utils.PROTEIN_DOMAIN
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

        #types of fusions: in-frame, out-of-frame, and any combination of
        # UTR, intron, intergenic, exonic (no known CDS),
        # CDS(not-reliable-start-or-end), CDS(truncated), CDS(complete)

        self.type=''

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
            try:
                if self.transcript_protein_junction_5prime < int(d[2]):
                    continue
                elif self.transcript_protein_junction_5prime >= int(d[3]):
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
            d[0]=self.name
            if self.transcript_protein_junction_3prime > int(d[3]):
                continue
            elif self.transcript_protein_junction_5prime <= int(d[2]):
                fusion_domains.append(d)
            else:
                d[2] = self.transcript_protein_junction_5prime
                fusion_domains.append(d)

        if self.name not in self.domains['pfam']:
            self.domains['pfam'][self.name] = [fusion_domains]
        else:
            self.domains['pfam'][self.name].append(fusion_domains)
        print self.domains['pfam'][self.name]

    def fetch_protein_seq(self):
        """
        Predict the potential protain amino acid sequence
        """

        self.protein_names = self.transcript1.protein_id + '-' + self.transcript2.protein_id

        self.transcript_protein_junction_5prime = int(self.transcript_cds_junction_5prime/3.)
        self.transcript_protein_junction_3prime = len(self.transcript2.protein_sequence) - int(self.transcript_cds_junction_3prime/3.)

        protein_seq = self.cds.seq.translate()
        protein_seq = protein_seq[0:protein_seq.find('*')]

        self.transcript_protein_molecular_weight = SeqUtils.molecular_weight(protein_seq)/1000.

        protein_seq = SeqRecord.SeqRecord(
            protein_seq,
            id=self.protein_names,
            name=self.protein_names,
            description="length=" + str(len(protein_seq)) + \
                ", kD: " + str(self.transcript_protein_molecular_weight) + \
                ", transcripts: " + str(self.name) + ', ' + \
                ", genes: " + str(self.gene_names)
        )

        self.protein = protein_seq

        self._annotate()

    def fetch_transcript_cds(self):
        """
        Predict the potential nucleotide sequence
        """

        transcript_seq = ''
        transcript_junction_5prime=0
        self.transcript_cds_junction_5prime=0
        self.transcript_cds_junction_3prime=0

        #5prime transcript

        if self.transcript1.strand=="+":
            for cds in self.transcript1.coding_sequence_position_ranges:
                if self.gene5prime.junction >= cds[1]:
                    self.transcript_cds_junction_5prime += cds[1] - cds[0]
                elif self.gene5prime.junction <= cds[0]:
                    break
                else:
                    self.transcript_cds_junction_5prime+=self.gene5prime.junction-cds[0]
                    break
        else:
            for cds in self.transcript1.coding_sequence_position_ranges:
                if self.gene5prime.junction <= cds[0]:
                    self.transcript_cds_junction_5prime += cds[1] - cds[0]
                elif self.gene5prime.junction >= cds[1]:
                    break
                else:
                    self.transcript_cds_junction_5prime += cds[1] - self.gene5prime.junction

        transcript_seq += self.transcript1.coding_sequence[0:self.transcript_cds_junction_5prime]

        if self.transcript2.strand=="+":
            for cds in self.transcript2.coding_sequence_position_ranges:
                if self.gene3prime.junction >= cds[1]:
                    self.transcript_cds_junction_3prime += cds[1] - cds[0]
                elif self.gene3prime.junction <= cds[0]:
                    break
                else:
                    self.transcript_cds_junction_3prime += self.gene3prime.junction - cds[0]
        else:
            for cds in transcript2.coding_sequence_position_ranges:
                if self.gene3prime.junction <= cds[0]:
                    self.transcript_cds_junction_3prime += cds[1] - cds[0]
                elif self.gene3prime.junction >= cds[1]:
                    break
                else:
                    self.transcript_cds_junction_3prime += cds[1] - self.gene3prime.junction

        transcript_seq += self.transcript1.coding_sequence[self.transcript_cds_junction_3prime::]

        transcript_seq = SeqRecord.SeqRecord(
            Seq.Seq(transcript_seq,generic_dna),
            id=self.name,
            name=self.name,
            description="length=" + str(len(transcript_seq)) + \
                ", genes: " + \
                str(self.transcript1.gene.name) + '-' + \
                str(self.transcript2.gene.name)
        )

        self.cds = transcript_seq

    def fetch_transcript_cdna(self):
        """
        Predict the potential nucleotide sequence
        """

        transcript_seq = ''
        self.transcript_cdna_junction_5prime=0
        self.transcript_cdna_junction_3prime=0

        #5prime transcript

        if self.transcript1.strand=="+":
            for exon in self.transcript1.exons:
                if self.gene5prime.junction >= exon.end:
                    self.transcript_cdna_junction_5prime += exon.end - exon.start
                elif self.gene5prime.junction <= exon.start:
                    break
                else:
                    self.transcript_cdna_junction_5prime+=self.gene5prime.junction-exon.start
                    break
        else:
            for exon in self.transcript1.exons:
                if self.gene5prime.junction <= exon.start:
                    self.transcript_cdna_junction_5prime += exon.end - exon.start
                elif self.gene5prime.junction >= exon.end:
                    break
                else:
                    self.transcript_cdna_junction_5prime += exon.end - self.gene5prime.junction

        transcript_seq += self.transcript1.sequence[0:self.transcript_cdna_junction_5prime]

        if self.transcript2.strand=="+":
            for exon in self.transcript2.exons:
                if self.gene3prime.junction >= exon.end:
                    self.transcript_cdna_junction_3prime += exon.end - exon.start
                elif self.gene3prime.junction <= exon.start:
                    break
                else:
                    self.transcript_cdna_junction_3prime += self.gene3prime.junction - exon.start
        else:
            for exon in self.transcript2.exons:
                if self.gene3prime.junction <= exon.start:
                    self.transcript_cdna_junction_3prime += exon.end - exon.start
                elif self.gene3prime.junction >= exon.end:
                    break
                else:
                    self.transcript_cdna_junction_3prime += exon.end - self.gene3prime.junction

        transcript_seq += self.transcript1.sequence[self.transcript_cdna_junction_3prime::]

        transcript_seq = SeqRecord.SeqRecord(
            Seq.Seq(transcript_seq,generic_dna),
            id=self.name,
            name=self.name,
            description="length=" + str(len(transcript_seq)) + \
                ", alternative transcript names: " + str(self.transcript1.name) + ', ' + \
                str(self.transcript2.name)
        )

        self.cdna = transcript_seq
