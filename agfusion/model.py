import itertools

from agfusion import utils, exceptions
import pandas
from PIL import *
import pyensembl
from Bio import Seq, SeqIO, SeqRecord
from Bio.Alphabet import generic_dna

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

        self.transcript_models = {}

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

        types_to_skip = ['retained_intron','processed_transcript']

        fout = open(out_file,'w')

        for combo in list(itertools.product(
                        self.gene5prime.gene.transcripts,
                        self.gene3prime.gene.transcripts)):
            transcript1 = combo[0]
            transcript2 = combo[1]

            if not transcript1.complete or not transcript2.complete:
                continue

            try:
                seq = self.fetch_transcript_cds(transcript1,transcript2)
            except:
                import pdb; pdb.set_trace()

            SeqIO.write(seq,fout,"fasta")

        fout.close()

    def save_protein_sequences(self,out_file):
        for combo in list(itertools.product(
                        self.gene5prime.gene.transcripts,
                        self.gene3prime.gene.transcripts)):
            transcript1 = combo[0]
            transcript2 = combo[1]

            r = self.fetch_protein_seq(t1,t2)

    def fetch_transcript_cdna(self,transcript1,transcript2):
        """
        Predict the potential nucleotide sequence
        """

        transcript_seq = ''
        transcript_junction_5prime=0
        transcript_junction_3prime=0

        #5prime transcript

        if transcript1.strand=="+":
            for exon in transcript1.exons:
                if self.gene5prime.junction >= exon.end:
                    transcript_junction_5prime += exon.end - exon.start
                elif self.gene5prime.junction <= exon.start:
                    break
                else:
                    transcript_junction_5prime+=self.gene5prime.junction-exon.start
                    break
        else:
            for exon in transcript1.exons:
                if self.gene5prime.junction <= exon.start:
                    transcript_junction_5prime += exon.end - exon.start
                elif self.gene5prime.junction >= exon.end:
                    break
                else:
                    transcript_junction_5prime += exon.end - self.gene5prime.junction

        transcript_seq += transcript1.sequence[0:transcript_junction_5prime]

        if transcript2.strand=="+":
            for exon in transcript2.exons:
                if self.gene3prime.junction >= exon.end:
                    transcript_junction_3prime += exon.end - exon.start
                elif self.gene3prime.junction <= exon.start:
                    break
                else:
                    transcript_junction_3prime += self.gene3prime.junction - exon.start
        else:
            for exon in transcript2.exons:
                if self.gene3prime.junction <= exon.start:
                    transcript_junction_3prime += exon.end - exon.start
                elif self.gene3prime.junction >= exon.end:
                    break
                else:
                    transcript_junction_3prime += exon.end - self.gene3prime.junction

        transcript_seq += transcript1.sequence[transcript_junction_3prime::]

        transcript_seq = SeqRecord.SeqRecord(
            Seq.Seq(transcript_seq,generic_dna),
            id=transcript1.id + '-' + transcript2.id,
            name=transcript1.id + '-' + transcript2.id,
            description="length=" + str(len(transcript_seq)) + \
                ", alternative transcript names: " + str(transcript1.name) + ', ' + \
                str(transcript2.name)
        )

        return transcript_seq

    def fetch_transcript_cds(self,transcript1,transcript2):
        """
        Predict the potential nucleotide sequence
        """

        transcript_seq = ''
        transcript_junction_5prime=0
        transcript_junction_3prime=0

        #5prime transcript

        if transcript1.strand=="+":
            for cds in transcript1.coding_sequence_position_ranges:
                if self.gene5prime.junction >= cds[1]:
                    transcript_junction_5prime += cds[1] - cds[0]
                elif self.gene5prime.junction <= cds[0]:
                    break
                else:
                    transcript_junction_5prime+=self.gene5prime.junction-cds[0]
                    break
        else:
            for cds in transcript1.coding_sequence_position_ranges:
                if self.gene5prime.junction <= cds[0]:
                    transcript_junction_5prime += cds[1] - cds[0]
                elif self.gene5prime.junction >= cds[1]:
                    break
                else:
                    transcript_junction_5prime += cds[1] - self.gene5prime.junction

        transcript_seq += transcript1.coding_sequence[0:transcript_junction_5prime]

        if transcript2.strand=="+":
            for cds in transcript2.coding_sequence_position_ranges:
                if self.gene3prime.junction >= cds[1]:
                    transcript_junction_3prime += cds[1] - cds[0]
                elif self.gene3prime.junction <= cds[0]:
                    break
                else:
                    transcript_junction_3prime += self.gene3prime.junction - cds[0]
        else:
            for cds in transcript2.coding_sequence_position_ranges:
                if self.gene3prime.junction <= cds[0]:
                    transcript_junction_3prime += cds[1] - cds[0]
                elif self.gene3prime.junction >= cds[1]:
                    break
                else:
                    transcript_junction_3prime += cds[1] - self.gene3prime.junction

        transcript_seq += transcript1.coding_sequence[transcript_junction_3prime::]

        transcript_seq = SeqRecord.SeqRecord(
            Seq.Seq(transcript_seq,generic_dna),
            id=transcript1.id + '-' + transcript2.id,
            name=transcript1.id + '-' + transcript2.id,
            description="length=" + str(len(transcript_seq)) + \
                ", alternative transcript names: " + str(transcript1.name) + ', ' + \
                str(transcript2.name)
        )

        return transcript_seq

    def fetch_protein_seq(self,transcript1,transcript2):
        """
        Predict the potential protain amino acid sequence
        """

        pass
