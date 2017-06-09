
class DataBaseError(Exception):
    def __init__(self,e):
        Exception.__init__(self,e)
        self.e = e

class GeneNotFound(Exception):
    def __init__(self,gene):
        Exception.__init__(self,"No entry in Ensembl found for {0}! Check its spelling and if you are using the right genome build.".format(gene))
        self.gene = gene

class GeneIDException3prime(Exception):
    def __init__(self,gene):
        Exception.__init__(self,"No Ensembl ID found for {0}! Check its spelling and if you are using the right genome build.".format(gene))
        self.gene = gene

class GeneIDException5prime(Exception):
    def __init__(self,gene):
        Exception.__init__(self,"No Ensembl ID found for {0}! Check its spelling and if you are using the right genome build.".format(gene))
        self.gene = gene

class TooManyGenesException(Exception):
    def __init__(self,gene,ids,build):
        Exception.__init__(self,"Multiple Ensembl IDs found matching {0}: {1} for genome {2}. Specify which Ensembl ID.".format(gene,', '.join(ids),build))
        self.gene = gene

class JunctionException5prime(Exception):
    def __init__(self):
        Exception.__init__(self,"Junction not within gene boundaries")

class JunctionException3prime(Exception):
    def __init__(self):
        Exception.__init__(self,"Junction not within gene boundaries")
