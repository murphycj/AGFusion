""" exceptions.py
Exceptions used for error handling.
"""


class DataBaseError(Exception):
    """_summary_

    Args:
        Exception (Exception): Base python exception.
    """

    def __init__(self, error):
        Exception.__init__(self, error)
        self.error = error


class GeneIDException(Exception):
    """Exception to handle case of no Ensemble gene ID for gene.

    Args:
        Exception (Exception): Base python exception.
    """

    def __init__(self, gene):
        Exception.__init__(
            self,
            f"No Ensembl ID found for {gene}! Check its spelling and if you are "
            "using the right genome build.",
        )


class TooManyGenesException(Exception):
    """Exception to handle case of too many genes.

    Args:
        Exception (Exception): Base python exception.
    """

    def __init__(self, gene, ids, build):
        ids = ", ".join(ids)
        Exception.__init__(
            self,
            f"Multiple Ensembl IDs found matching {gene}: {ids} for genome {build}."
            " Specify which Ensembl ID.",
        )


class JunctionException(Exception):
    """Exception to handle junction not being in the gene.

    Args:
        Exception (Exception): Base python exception.
    """

    def __init__(self, gene, junction):
        Exception.__init__(self, f"Junction {junction} not within {gene} gene boundaries!")
