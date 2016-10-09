
class DataBaseError(Exception):
    def __init__(self,e):
        Exception.__init__(self,e)
        self.e = e

class GeneIDException(Exception):
    def __init__(self,gene):
        Exception.__init__(self,"No Ensembl ID found for {0}! Check to make sure you are using the right genome build.".format(gene))
        self.gene = gene

class TooManyGenesException(Exception):
    def __init__(self,gene):
        Exception.__init__(self,"Multiple ensembl IDs found matching {0}. Try specifying just the Ensembl ID.".format(gene))
        self.gene = gene

class JunctionException(Exception):
    def __init__(self):
        Exception.__init__(self,"Junction not within gene boundaries")
