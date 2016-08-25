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

class Protein:
    """
    Stores the necessary information for proteins
    """

    def __init__(self):
        self.id=''
        self.amino_acid_seq=''
        self.domains={
            i:{j:0,k:0} for i,j,k in utils.PROTEIN_DOMAIN
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

class Fusion(Model):
    """
    Generates the information needed for the gene fusion
    """

    def __init__(self,gene5prime=None,gene3prime=None,db=None):
        self.gene5prime=gene5prime
        self.gene3prime=gene3prime
        self.db=db

        self.transcript_cds = {}
        self.transcript_cdna = {}
        self.transcript_protein = {}
        self.transcript_protein_molecular_weight = {}

        self.transcript_cdna_junction_5prime = {}
        self.transcript_cdna_junction_3prime = {}
        self.transcript_cds_junction_5prime = {}
        self.transcript_cds_junction_3prime = {}
        self.transcript_protein_junction_5prime = {}
        self.transcript_protein_junction_3prime = {}

    def _annotate(self,t1,t2):
        """
        Annotate the gene fusion's protein using the protein annotaiton
        from its two genes
        """
        self.db.c.execute(
            'SELECT * FROM pfam WHERE ensembl_transcript_id==\"' + \
            t1.id + \
            '\"'
        )

        print(self.db.c.fetchall())

        self.db.c.execute(
            'SELECT * FROM pfam WHERE ensembl_transcript_id==\"' + \
            t1.id + \
            '\"'
        )
        print(self.db.c.fetchall())

    def save_transcript_cdna(self,out_file):

        fout = open(out_file,'w')

        for combo in list(itertools.product(
                        self.gene5prime.gene.transcripts,
                        self.gene3prime.gene.transcripts)):
            transcript1 = combo[0]
            transcript2 = combo[1]

            seq = self.fetch_transcript_cdna(transcript1,transcript2)

            SeqIO.write(seq,fout,"fasta")

        fout.close()

    def save_transcript_cds(self,out_file):

        fout = open(out_file,'w')

        for combo in list(itertools.product(
                        self.gene5prime.gene.transcripts,
                        self.gene3prime.gene.transcripts)):
            transcript1 = combo[0]
            transcript2 = combo[1]

            if not transcript1.complete or not transcript2.complete:
                continue

            seq = self.fetch_transcript_cds(transcript1,transcript2)

            SeqIO.write(seq,fout,"fasta")

        fout.close()

    def save_proteins(self,out_file):

        fout = open(out_file,'w')

        for combo in list(itertools.product(
                        self.gene5prime.gene.transcripts,
                        self.gene3prime.gene.transcripts)):
            transcript1 = combo[0]
            transcript2 = combo[1]

            if not transcript1.complete or not transcript2.complete:
                continue

            name = transcript1.id + '-' + transcript2.id

            if name not in self.transcript_cds:
                self.fetch_transcript_cds(transcript1,transcript2)

            seq = self.fetch_protein_seq(transcript1,transcript2)

            SeqIO.write(seq,fout,"fasta")

        fout.close()

    def fetch_transcript_cdna(self,transcript1,transcript2):
        """
        Predict the potential nucleotide sequence
        """

        transcript_seq = ''
        self.transcript_cdna_junction_5prime=0
        self.transcript_cdna_junction_3prime=0

        #5prime transcript

        if transcript1.strand=="+":
            for exon in transcript1.exons:
                if self.gene5prime.junction >= exon.end:
                    self.transcript_cdna_junction_5prime += exon.end - exon.start
                elif self.gene5prime.junction <= exon.start:
                    break
                else:
                    self.transcript_cdna_junction_5prime+=self.gene5prime.junction-exon.start
                    break
        else:
            for exon in transcript1.exons:
                if self.gene5prime.junction <= exon.start:
                    self.transcript_cdna_junction_5prime += exon.end - exon.start
                elif self.gene5prime.junction >= exon.end:
                    break
                else:
                    self.transcript_cdna_junction_5prime += exon.end - self.gene5prime.junction

        transcript_seq += transcript1.sequence[0:self.transcript_cdna_junction_5prime]

        if transcript2.strand=="+":
            for exon in transcript2.exons:
                if self.gene3prime.junction >= exon.end:
                    self.transcript_cdna_junction_3prime += exon.end - exon.start
                elif self.gene3prime.junction <= exon.start:
                    break
                else:
                    self.transcript_cdna_junction_3prime += self.gene3prime.junction - exon.start
        else:
            for exon in transcript2.exons:
                if self.gene3prime.junction <= exon.start:
                    self.transcript_cdna_junction_3prime += exon.end - exon.start
                elif self.gene3prime.junction >= exon.end:
                    break
                else:
                    self.transcript_cdna_junction_3prime += exon.end - self.gene3prime.junction

        transcript_seq += transcript1.sequence[self.transcript_cdna_junction_3prime::]

        name = transcript1.id + '-' + transcript2.id

        transcript_seq = SeqRecord.SeqRecord(
            Seq.Seq(transcript_seq,generic_dna),
            id=name,
            name=name,
            description="length=" + str(len(transcript_seq)) + \
                ", alternative transcript names: " + str(transcript1.name) + ', ' + \
                str(transcript2.name)
        )

        self.transcript_cdna[name] = transcript_seq

        return transcript_seq

    def fetch_transcript_cds(self,transcript1,transcript2):
        """
        Predict the potential nucleotide sequence
        """

        transcript_seq = ''
        transcript_junction_5prime=0
        self.transcript_cds_junction_5prime=0
        self.transcript_cds_junction_3prime=0

        #5prime transcript

        if transcript1.strand=="+":
            for cds in transcript1.coding_sequence_position_ranges:
                if self.gene5prime.junction >= cds[1]:
                    self.transcript_cds_junction_5prime += cds[1] - cds[0]
                elif self.gene5prime.junction <= cds[0]:
                    break
                else:
                    self.transcript_cds_junction_5prime+=self.gene5prime.junction-cds[0]
                    break
        else:
            for cds in transcript1.coding_sequence_position_ranges:
                if self.gene5prime.junction <= cds[0]:
                    self.transcript_cds_junction_5prime += cds[1] - cds[0]
                elif self.gene5prime.junction >= cds[1]:
                    break
                else:
                    self.transcript_cds_junction_5prime += cds[1] - self.gene5prime.junction

        transcript_seq += transcript1.coding_sequence[0:self.transcript_cds_junction_5prime]

        if transcript2.strand=="+":
            for cds in transcript2.coding_sequence_position_ranges:
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

        transcript_seq += transcript1.coding_sequence[self.transcript_cds_junction_3prime::]

        name = transcript1.id + '-' + transcript2.id

        transcript_seq = SeqRecord.SeqRecord(
            Seq.Seq(transcript_seq,generic_dna),
            id=transcript1.id + '-' + transcript2.id,
            name=transcript1.id + '-' + transcript2.id,
            description="length=" + str(len(transcript_seq)) + \
                ", genes: " + \
                str(transcript1.gene.name) + '-' + \
                str(transcript2.gene.name)
        )

        self.transcript_cds[name] = transcript_seq

        return transcript_seq

    def fetch_protein_seq(self,transcript1,transcript2):
        """
        Predict the potential protain amino acid sequence
        """

        protein_names = transcript1.protein_id + '-' + transcript2.protein_id
        transcript_names = transcript1.id + '-' + transcript2.id
        gene_names = transcript1.gene.name + '-' + transcript2.gene.name


        protein_seq = \
            transcript1.protein_sequence[
                0:self.transcript_cds_junction_5prime
            ] + \
            transcript2.protein_sequence[
                (len(transcript2.protein_sequence)-self.transcript_cds_junction_3prime)::
            ]

        protein_seq = Seq.Seq(protein_seq,generic_protein)
        self.transcript_protein_molecular_weight[transcript_names] = SeqUtils.molecular_weight(protein_seq)/1000.

        protein_seq = SeqRecord.SeqRecord(
            protein_seq,
            id=protein_names,
            name=protein_names,
            description="length=" + str(len(protein_seq)) + \
                ", kD: " + str(self.transcript_protein_molecular_weight[transcript_names]) + \
                ", transcripts: " + str(transcript_names) + ', ' + \
                ", genes: " + str(gene_names)
        )

        self.transcript_protein[transcript_names] = protein_seq

        return protein_seq
