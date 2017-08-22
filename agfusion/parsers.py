"""
Parses output files from fusion-finding algorithms
"""

import re

class _Parser(object):
    def __init__(self):
        self.fusions = []
        self.iterator = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterator >= len(self.fusions):
            raise StopIteration
        else:
            self.iterator += 1
            return self.fusions[self.iterator-1]

    next = __next__

class STARFusion(_Parser):
    def __init__(self,infile):
        super(STARFusion, self).__init__()

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

class FusionEntry():
    pass

class Bellerophontes(_Parser):
    pass

class Chimerascan(_Parser):
    pass

class EricScript(_Parser):
    pass

class FusionCatcher(_Parser):
    def __init__(self,infile):
        super(FusionCatcher, self).__init__()

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

class FusionHunter(_Parser):
    def __init__(self,infile):
        super(FusionHunter, self).__init__()

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


class FusionMap(_Parser):
    def __init__(self,infile):
        super(FusionMap, self).__init__()

        fin = open(infile,'r')
        for line in fin.readlines():
            if re.findall('^FusionID',line):
                line = line.strip().split('\t')
                for i, j in zip([6,8,9,13],['Position1','Position2','KnownGene1','KnownGene2']):
                    assert line[i]==j, "Unrecognized FusionMap input: %s".format(j)
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

class JAFFA(_Parser):
    pass

class MapSplice(_Parser):
    pass

class nFuse(_Parser):
    pass

class SOAPfuse(_Parser):
    pass

class TopHatFusion(_Parser):
    def __init__(self,infile):
        super(TopHatFusion, self).__init__()

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

parsers = {
    #'bellerophontes':Bellerophontes,
    #'chimerascan':Chimerascan,
    #'ericscript':EricScript,
    'fusioncatcher':FusionCatcher,
    'fusionhunter':FusionHunter,
    'fusionmap':FusionMap,
    #'jaffa':JAFFA,
    #'mapsplice':MapSplice,
    #'nfuse':nFuse,
    #'soapfuse':SOAPfuse,
    'tophatfusion':TopHatFusion,
    'starfusion':STARFusion
}
