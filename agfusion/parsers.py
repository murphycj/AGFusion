"""
Parses output files from fusion-finding algorithms
"""

import re

class _Parser(object):
    def __init__(self, logger):
        self.fusions = []
        self.iterator = 0
        self.logger = logger

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterator >= len(self.fusions):
            raise StopIteration
        else:
            self.iterator += 1
            return self.fusions[self.iterator-1]

    def _check_data(self):
        if len(self.fusions)==0:
            self.logger.error("Read 0 fusions from the file! Exiting...")
        else:
            self.logger.info(
                "Read {} fusions from the file.".format(len(self.fusions))
            )
    next = __next__

class STARFusion(_Parser):
    def __init__(self,infile,logger):
        super(STARFusion, self).__init__(logger)

        fin = open(infile,'r')
        for line in fin.readlines():
            if re.findall('^#',line):
                line = line.rstrip().split('\t')
                assert line[0]=='#FusionName', 'Unrecognized STAR-Fusion input'
                assert line[4]=='LeftGene', 'Unrecognized STAR-Fusion input'
                assert line[5]=='LeftBreakpoint', 'Unrecognized STAR-Fusion input'
                assert line[6]=='RightGene', 'Unrecognized STAR-Fusion input'
                assert line[7]=='RightBreakpoint', 'Unrecognized STAR-Fusion input'
                continue

            line = line.strip().split('\t')

            gene_5prime = line[4].split('^')[1].split('.')[0]
            gene_5prime_name = line[4].split('^')[0]
            gene_5prime_junction = int(line[5].split(':')[1])
            gene_3prime = line[6].split('^')[1].split('.')[0]
            gene_3prime_name = line[6].split('^')[0]
            gene_3prime_junction = int(line[7].split(':')[1])
            self.fusions.append(
                {
                    'ensembl_5prime':gene_5prime,
                    'ensembl_3prime':gene_3prime,
                    'alternative_name_5prime':gene_5prime_name,
                    'alternative_name_3prime':gene_3prime_name,
                    'junction_5prime':gene_5prime_junction,
                    'junction_3prime':gene_3prime_junction
                }
            )
        fin.close()

        self._check_data()

class EricScript(_Parser):
    def __init__(self,infile,logger):
        super(EricScript, self).__init__(logger)

        fin = open(infile,'r')
        for line in fin.readlines():
            if re.findall('^GeneName1',line):
                line = line.strip().split('\t')
                for i, j in zip([3,6,8,9],['Breakpoint1','Breakpoint2','EnsemblGene1','EnsemblGene2']):
                    assert line[i]==j, "Unrecognized EricScript input: {}".format(j)
            else:
                line = line.strip().split('\t')

                gene_5prime_name = line[8]
                gene_5prime_junction = int(line[3])
                gene_3prime_name = line[9]
                gene_3prime_junction = int(line[6])
                self.fusions.append(
                    {
                        'ensembl_5prime':None,
                        'ensembl_3prime':None,
                        'alternative_name_5prime':gene_5prime_name,
                        'alternative_name_3prime':gene_3prime_name,
                        'junction_5prime':gene_5prime_junction,
                        'junction_3prime':gene_3prime_junction
                    }
                )
        fin.close()

        self._check_data()

class FusionCatcher(_Parser):
    def __init__(self,infile,logger):
        super(FusionCatcher, self).__init__(logger)

        fin = open(infile,'r')
        for line in fin.readlines():
            if re.findall('^Gene_1_symbol',line):
                line = line.rstrip().split('\t')
                assert line[8]=='Fusion_point_for_gene_1(5end_fusion_partner)', 'Unrecognized FusionCatcher input'
                assert line[9]=='Fusion_point_for_gene_2(3end_fusion_partner)', 'Unrecognized FusionCatcher input'
                assert line[10]=='Gene_1_id(5end_fusion_partner)', 'Unrecognized FusionCatcher input'
                assert line[11]=='Gene_2_id(3end_fusion_partner)', 'Unrecognized FusionCatcher input'
                continue

            line = line.strip().split('\t')
            self.fusions.append(
                {
                    'ensembl_5prime':line[10],
                    'ensembl_3prime':line[11],
                    'alternative_name_5prime':line[0],
                    'alternative_name_3prime':line[1],
                    'junction_5prime':int(line[8].split(':')[1]),
                    'junction_3prime':int(line[9].split(':')[1])
                }
            )
        fin.close()

        self._check_data()

class FusionHunter(_Parser):
    def __init__(self,infile,logger):
        super(FusionHunter, self).__init__(logger)

        gene1 = gene2 = None
        gene1_junction = gene2_junction = None

        fin = open(infile,'r')
        for line in fin.readlines():

            if re.findall('^# Fusion:',line):

                if gene1 is not None and gene2 is not None:
                    self.fusions.append(
                        {
                            'ensembl_5prime':None,
                            'ensembl_3prime':None,
                            'alternative_name_5prime':gene1,
                            'alternative_name_3prime':gene2,
                            'junction_5prime':int(gene1_junction),
                            'junction_3prime':int(gene2_junction)
                        }
                    )

                strands = re.findall("(?<=\[)(.*)(?=\])",line)
                assert len(strands)==1,"Unrecognized FusionHunter input. Incorrect strand information."
                gene1_strand = strands[0][0]
                gene2_strand = strands[0][1]

            elif re.findall('^--',line):
                #new breakpoint
                if gene1 is not None and gene2 is not None:
                    self.fusions.append(
                        {
                            'ensembl_5prime':None,
                            'ensembl_3prime':None,
                            'alternative_name_5prime':gene1,
                            'alternative_name_3prime':gene2,
                            'junction_5prime':int(gene1_junction),
                            'junction_3prime':int(gene2_junction)
                        }
                    )

            elif re.findall('^->',line):
                junctions = re.findall("(chr[0-9]*):([0-9]*)-([0-9]*)",line)
                assert len(junctions)==2,"Unrecognized FusionHunter input. Incorrect junction information."

                gene1 = gene2 = None
                gene1_junction = gene2_junction = None

                gene1, gene2 = re.findall("[A-Z0-9]*\sx\s[A-Z0-9a-z]*(?=\t)",line)[0].split(' x ')
                assert gene1 is not None and gene2 is not None, "Unrecognized FusionHunter input. Incorrect gnee information."

                if gene1_strand=='+':
                    gene1_junction = junctions[0][2]
                else:
                    gene1_junction = junctions[0][1]

                if gene2_strand=='+':
                    gene2_junction = junctions[1][1]
                else:
                    gene2_junction = junctions[1][2]
        fin.close()

        if gene1 is not None and gene2 is not None:
            self.fusions.append(
                {
                    'ensembl_5prime':None,
                    'ensembl_3prime':None,
                    'alternative_name_5prime':gene1,
                    'alternative_name_3prime':gene2,
                    'junction_5prime':int(gene1_junction),
                    'junction_3prime':int(gene2_junction)
                }
            )

        self._check_data()


class FusionMap(_Parser):
    def __init__(self,infile,logger):
        super(FusionMap, self).__init__(logger)

        fin = open(infile,'r')
        for line in fin.readlines():
            if re.findall('^FusionID',line):
                line = line.strip().split('\t')
                for i, j in zip([6,8,9,13],['Position1','Position2','KnownGene1','KnownGene2']):
                    assert line[i]==j, "Unrecognized FusionMap input {} header not in expected position.".format(j)
            else:
                line = line.strip().split('\t')

                gene_5prime_name = line[9]
                gene_5prime_junction = int(line[6])
                gene_3prime_name = line[13]
                gene_3prime_junction = int(line[8])
                self.fusions.append(
                    {
                        'ensembl_5prime':None,
                        'ensembl_3prime':None,
                        'alternative_name_5prime':gene_5prime_name,
                        'alternative_name_3prime':gene_3prime_name,
                        'junction_5prime':gene_5prime_junction,
                        'junction_3prime':gene_3prime_junction
                    }
                )
        fin.close()

        self._check_data()

class MapSplice(_Parser):
    def __init__(self,infile,logger):
        super(MapSplice, self).__init__(logger)

        fin = open(infile,'r')
        for line in fin.readlines():
            if re.findall('^chrom',line):
                line = line.strip().split('\t')
                for i, j in zip([1,2,60,61],['doner_end','acceptor_start','annotated_gene_donor','annotated_gene_acceptor']):
                    assert line[i]==j, "Unrecognized MapSplice input. {} header not in expected position.".format(j)
            else:
                line = line.strip().split('\t')

                gene_5prime_name = line[60]
                gene_5prime_junction = int(line[1])
                gene_3prime_name = line[61]
                gene_3prime_junction = int(line[2])
                self.fusions.append(
                    {
                        'ensembl_5prime':None,
                        'ensembl_3prime':None,
                        'alternative_name_5prime':gene_5prime_name,
                        'alternative_name_3prime':gene_3prime_name,
                        'junction_5prime':gene_5prime_junction,
                        'junction_3prime':gene_3prime_junction
                    }
                )
        fin.close()

        self._check_data()

class TopHatFusion(_Parser):
    def __init__(self,infile,logger):
        super(TopHatFusion, self).__init__(logger)

        fin = open(infile,'r')
        for line in fin.readlines():
            line = line.strip().split('\t')

            gene_5prime = line[1]
            gene_5prime_name = line[1]
            gene_5prime_junction = int(line[3])
            gene_3prime = line[4]
            gene_3prime_name = line[4]
            gene_3prime_junction = int(line[6])
            self.fusions.append(
                {
                    'ensembl_5prime':gene_5prime,
                    'ensembl_3prime':gene_3prime,
                    'alternative_name_5prime':gene_5prime_name,
                    'alternative_name_3prime':gene_3prime_name,
                    'junction_5prime':gene_5prime_junction,
                    'junction_3prime':gene_3prime_junction
                }
            )
        fin.close()

        self._check_data()

class DeFuse(_Parser):
    def __init__(self,infile,logger):
        super(DeFuse, self).__init__(logger)

        data_indices = {
            'gene5prime':None,
            'gene3prime':None,
            'gene5prime_junction':None,
            'gene3prime_junction':None
        }

        fin = open(infile,'r')
        for line in fin.readlines():
            if re.findall('^cluster_id',line):
                line = line.strip().split('\t')
                for column in data_indices.keys():
                    try:
                        data_indices[column] = line.index(column)
                    except ValueError:
                        logger.error(
                            "Unrecognized {} input! Cannot find {} column."
                            .format(
                                self.__class__.__name__,
                                column
                            )
                        )
                        exit()
                continue
            if not any([i is None for i in data_indices.keys()]):
                line = line.strip().split('\t')
                self.fusions.append(
                    {
                        'gene5prime':line[data_indices['gene5prime']],
                        'gene3prime':line[data_indices['gene3prime']],
                        'gene5prime_junction':int(line[data_indices['gene5prime_junction']]),
                        'gene3prime_junction':int(line[data_indices['gene3prime_junction']])
                    }
                )
        fin.close()

        self._check_data()

class Chimerascan(_Parser):
    def __init__(self,infile,logger):
        super(Chimerascan, self).__init__(logger)

        data_indices = {
            'gene5prime':None,
            'gene3prime':None,
            'gene5prime_junction':None,
            'gene3prime_junction':None,
            'gene5prime_strand':None,
            'gene3prime_strand':None
        }

        fin = open(infile,'r')
        for line in fin.readlines():
            if re.findall('^#chrom5p',line):
                line = line.strip().split('\t')
                for column in data_indices.keys():
                    try:
                        data_indices[column] = line.index(column)
                    except ValueError:
                        logger.error(
                            "Unrecognized {} input! Cannot find {} column."
                            .format(
                                self.__class__.__name__,
                                column
                            )
                        )
                        exit()
                continue
            elif gene1_ix is not None:
                line = line.strip().split('\t')

                gene_5prime_name = line[gene1_ix]
                gene_5prime_junction = int(line[gene1_junction_ix])
                gene_3prime_name = line[gene2_ix]
                gene_3prime_junction = int(line[gene2_junction_ix])
                self.fusions.append(
                    {
                        'gene5prime':gene_5prime_name,
                        'gene3prime':gene_3prime_name,
                        'gene5prime_junction':gene_5prime_junction,
                        'gene3prime_junction':gene_3prime_junction
                    }
                )
        fin.close()

        self._check_data()

class ChimeRScope(_Parser):
    def __init__(self,infile,logger):
        super(ChimeRScope, self).__init__(logger)

        data_indices = {
            'Gene1':2,
            'Gene2':4,
            'Gene1_fusionPoint':7,
            'Gene2_fusionPoint':9
        }

        fin = open(infile,'r')
        for line in fin.readlines():
            if re.findall('^ConfidentScore',line):
                line = line.strip().split('\t')
                for column in data_indices.keys():
                    try:
                        data_indices[column] = line.index(column)
                    except ValueError:
                        logger.error(
                            "Unrecognized {} input! Cannot find {} column."
                            .format(
                                self.__class__.__name__,
                                column
                            )
                        )
                        exit()
                continue

            line = line.strip().split('\t')
            self.fusions.append(
                {
                    'gene5prime':line[data_indices['Gene1']],
                    'gene3prime':line[data_indices['Gene2']],
                    'gene5prime_junction':line[data_indices['Gene1_fusionPoint']],
                    'gene3prime_junction':line[data_indices['Gene2_fusionPoint']]
                }
            )
        fin.close()

        self._check_data()

class JAFFA(_Parser):
    def __init__(self,infile,logger):
        super(JAFFA, self).__init__(logger)

        data_indices = {
            'base1':7,
            'base2':9,
            'fusion genes':1
        }

        fin = open(infile,'r')
        for line in fin.readlines():
            if re.findall('sample',line):
                line = line.strip().replace('"','').split(',')
                for column in data_indices.keys():
                    try:
                        data_indices[column] = line.index(column)
                    except ValueError:
                        logger.error(
                            "Unrecognized {} input! Cannot find {} column."
                            .format(
                                self.__class__.__name__,
                                column
                            )
                        )
                        exit()
                continue

            line = line.strip().replace('"','').split(',')
            self.fusions.append(
                {
                    'gene5prime':line[data_indices['fusion genes']].split(':')[0],
                    'gene3prime':line[data_indices['fusion genes']].split(':')[1],
                    'gene5prime_junction':int(line[data_indices['base1']]),
                    'gene3prime_junction':int(line[data_indices['base2']])
                }
            )
        fin.close()

        self._check_data()

class Bellerophontes(_Parser):
    def __init__(self,infile,logger):
        super(Bellerophontes, self).__init__(logger)

        fin = open(infile,'r')
        for line in fin.readlines():

            line = line.strip().split('\t')
            if len(line) <= 2:
                continue

            gene1_strand = line[2]
            gene2_strand = line[6]

            self.fusions.append(
                {
                    'gene5prime':line[0],
                    'gene3prime':line[4],
                    'gene5prime_junction':int(line[9]),
                    'gene3prime_junction':int(line[11])
                }
            )
        fin.close()

        self._check_data()

class NFuse(_Parser):
    def __init__(self,infile,logger):
        super(NFuse, self).__init__(logger)

        self._check_data()

class SOAPfuse(_Parser):
    def __init__(self,infile,logger):
        super(SOAPfuse, self).__init__(logger)

        self._check_data()

class FusionEntry():
    def __init__(self,infile,logger):
        super(FusionEntry, self).__init__(logger)

        self._check_data()

parsers = {
    'bellerophontes':Bellerophontes,
    #'chimerascan':Chimerascan,
    'chimerscope':ChimeRScope,
    'defuse':DeFuse,
    'ericscript':EricScript,
    'fusioncatcher':FusionCatcher,
    'fusionhunter':FusionHunter,
    'fusionmap':FusionMap,
    'jaffa':JAFFA,
    #'fusionseq':FusionSeq,
    #'prada':Prada,
    #'gfusion':GFusion,
    #'machete':Machete,
    'mapsplice':MapSplice,
    #'nfuse':nFuse,
    #'soapfuse':SOAPfuse,
    'starfusion':STARFusion,
    'tophatfusion':TopHatFusion
}
