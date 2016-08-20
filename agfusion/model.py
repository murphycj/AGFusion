from agfusion import utils
from PIL import *

class Model(object):
    def __new__(cls, *args, **kwargs):
        if cls is Model:
            raise TypeError("Model shoud not be instantiated")
        return object.__new__(cls, *args, **kwargs)

    def save_image(self,):
        pass


class Gene(Model):
    """
    Stores the necessary information to specify the architecture of either
    wild-type genes or fusion gene.
    """

    def __init__(self,ensembl_gene_id='',junction=0,db=None):
        self.ensembl_gene_id=ensembl_gene_id
        self.junction=junction
        self.symbol=''
        self.ensembl_transcript_ids=[]
        self.chr=''
        self.start=0
        self.end=0
        self.db=db

        self._fetch()

    def _fetch(self):
        """
        Fetch the appropriate information from the SQLite3 database
        """

        pass

class Transcript:
    """
    Stores the necessary information for transcript
    """

    def __init__(self):
        self.id=''
        self.biotype=''
        self.seq=''
        self.length=0
        self.ensembl_protein_id=''

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

    def __init__(self,gene5prime=None,gene3prime=None):
        self.gene5prime=gene5prime
        self.gene3prime=gene3prime

        self.junction_protein=0

        self.transcripts=[]

    def _create_transcripts(self):
        pass

    def _create_proteins(self):
        pass
