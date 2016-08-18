
class Gene:
    """
    Stores the necessary information to specify the architecture of either
    wild-type genes or fusion gene.
    """

    def __init__(self):
        self.gene_5prime=''
        self.gene_5prime_junction=0

        self.gene_3prime=''
        self.gene_3prime_junction=0
