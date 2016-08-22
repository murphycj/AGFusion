from agfusion import utils, exceptions
import pandas
from PIL import *
import pyfaidx

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

    def __init__(self,query_gene='',junction=0,reference='',db=None):
        self.query_gene = query_gene
        self.reference=reference
        self.db=db
        self.junction=junction

        self._fetch_and_validate()

        self.symbol=''
        self.entrez=0
        self.ensembl_transcript={}
        self.chr=''
        self.start=0
        self.end=0

        self._fetch_transcripts()

    def _fetch_chr(self):
        pass

    def _fetch_and_validate(self):
        """
        Fetch chromosome, strand information. Also confirms the query gene
        name is either ensembl gene id, gene symbol, or ensemble gene id.
        Also validate that the fusion junction is within the gene boundary
        """

        self.ensembl_gene_id = self.query_gene

        self.db.c.execute(
            "SELECT * FROM " + self.reference + \
            " WHERE ensembl_gene_id==\"" + self.ensembl_gene_id + "\""
        )

        d = pandas.DataFrame(self.db.c.fetchall(),columns=utils.GENE_SCHEMA_COLUMNS)

        assert d.shape[0]==1, exceptions.DataBaseError(
            "Too many or too few ensembl_gene_id records were" + \
            " retrieved for %s" % self.ensembl_gene_id
        )

        self.chr = d.ix[0]['chr']
        self.strand = d.ix[0]['strand']

        if self.junction > d.ix[0]['genomic_end'] or self.junction < d.ix[0]['genomic_start']:
            raise exceptions.JunctionException()

    def _fetch_transcripts(self):
        """
        Fetches the ensembl transcript ids as well as exon and CDS coordinates
        from SQLite3 database
        """

        #get transcript information for the gene

        self.db.c.execute(
            "SELECT * FROM " + \
            self.reference + "_transcript WHERE " + \
            "ensembl_gene_id == \"" + self.ensembl_gene_id + "\""
        )

        transcripts = pandas.DataFrame(self.db.c.fetchall(),columns=utils.TRANSCRIPT_SCHEMA_COLUMNS)

        assert transcripts.shape[0]>0, exceptions.DataBaseError(
            "Could not find transcript information for" + \
            " %s" % self.ensembl_gene_id
        )

        for i in transcripts.index:

            #get annotaiton information for each transcript

            self.db.c.execute(
                "SELECT * FROM " + \
                self.reference + "_annotation_transcript WHERE " + \
                "ensembl_transcript_id == \"" + \
                transcripts.ix[i]['ensembl_transcript_id'] + \
                "\""
            )

            transcript_annotation = pandas.DataFrame(
                self.db.c.fetchall(),
                columns=utils.TRANSCRIPT_ANNOTATION_COLUMNS
            )

            assert transcript_annotation.shape[0]>0, exceptions.DataBaseError(
                "Could not find transcript annotation information for" + \
                " %s" % transcripts.ix[i]['ensembl_transcript_id']
            )

            transcript = Transcript(
                ensembl_transcript_id=transcripts.ix[i]['ensembl_transcript_id'],
                ensembl_protein_id=transcripts.ix[i]['ensembl_protein_id'],
                biotype=transcripts.ix[i]['transcript_biotype'],
                start=transcripts.ix[i]['transcript_genomic_start'],
                end=transcripts.ix[i]['transcript_genomic_end'],
                annotation=transcript_annotation
            )

            self.ensembl_transcript[transcripts.ix[i]['ensembl_transcript_id']] = transcript


    def predict_effect(self):
        """
        By examining what feature (e.g. exon, CDS, UTR, ...) the fusion junction
        is located for a give transcript determines the predicted effect
        on the gene
        """

        for ensembl_transcript_id, transcript in self.ensembl_transcript.items():
            transcript.predict_effect(self.junction)


class Transcript:
    """
    Stores the necessary information for transcript
    """

    def __init__(self,ensembl_transcript_id,ensembl_protein_id,biotype,start,end,annotation):
        self.ensembl_transcript_id=ensembl_transcript_id
        self.biotype=biotype
        self.ensembl_protein_id=ensembl_protein_id
        self.genomic_start=start
        self.genomic_end=end
        self.annotation=annotation
        self.effect=''

    def predict_effect(self,junction):
        """
        By examining what feature (e.g. exon, CDS, UTR, ...) the fusion junction
        is located within the transcript determines the predicted effect
        """
        import pdb; pdb.set_trace()


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

    def _create_transcript_seq(self):
        """
        Predict the potential nucleotide sequence
        """

        pass

    def _create_protein_seq(self):
        """
        Predict the potential protain amino acid sequence
        """

        pass

    def _annotate(self):
        """
        Annotate the gene fusion's protein using the protein annotaiton
        from its two genes
        """
        pass
