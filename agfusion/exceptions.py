
class DataBaseError(Exception):
    def __init__(self,e):
        Exception.__init__(self,e)
        self.e = e

class GeneIDException(Exception):
    def __init__(self,gene):
        Exception.__init__(
            self,
            "No Ensembl ID found for {}! Check its spelling and if you are " \
            "using the right genome build.".format(gene)
        )

class TooManyGenesException(Exception):
    def __init__(self,gene,ids,build):
        Exception.__init__(
            self,
            "Multiple Ensembl IDs found matching {}: {} for genome {}." \
            " Specify which Ensembl ID.".format(gene,', '.join(ids),build)
        )

class JunctionException(Exception):
    def __init__(self,gene,junction):
        Exception.__init__(
            self,
            "Junction {} not within {} gene boundaries!"
            .format(junction,gene)
        )
