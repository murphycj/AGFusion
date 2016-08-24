import itertools

from agfusion import utils, exceptions
import pandas
from PIL import *
import pyensembl
from Bio import Seq

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
        self.db.c.execute('SELECT * FROM pfam WHERE ensembl_transcript_id==\"' + t1.id + '\"')

        print(self.db.c.fetchall())

        self.db.c.execute('SELECT * FROM pfam WHERE ensembl_transcript_id==\"' + t1.id + '\"')
        print(self.db.c.fetchall())

    def save_transcript_sequences(self,out_file):
        for combo in list(itertools.product(self.gene5prime.gene.transcripts,self.gene3prime.gene.transcripts)):
            t1 = combo[0]
            t2 = combo[1]

            r = self.fetch_transcript_seq(t1,t2)

    def save_protein_sequences(self,out_file):
        for combo in list(itertools.product(self.gene5prime.gene.transcripts,self.gene3prime.gene.transcripts)):
            t1 = combo[0]
            t2 = combo[1]

            r = self.fetch_protein_seq(t1,t2)

    def fetch_transcript_cdna(self,transcript1,transcript2):
        """
        Predict the potential nucleotide sequence
        """

        cdna_seq = ''

        #5prime transcript

        length_5prime=0

        if transcript1.strand=="+":
            for exon in transcript1.exons:
                if self.gene5prime.junction >= exon.end:
                    length_5prime+=exon.end - exon.end
                elif self.gene5prime.junction <= exon.start:
                    break
                else:
                    length_5prime+=self.gene5prime.junction-exon.start
                    break
        else:


            


    def fetch_transcript_cds(self,transcript1,transcript2):
        """
        Predict the potential nucleotide sequence
        """

        pass

    def fetch_protein_seq(self,transcript1,transcript2):
        """
        Predict the potential protain amino acid sequence
        """

        pass
