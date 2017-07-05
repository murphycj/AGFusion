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
    pass

class FusionMap(_Parser):
    pass

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
    'bellerophontes':Bellerophontes,
    'chimerascan':Chimerascan,
    'ericscript':EricScript,
    'fusioncatcher':FusionCatcher,
    'fusionhunter':FusionHunter,
    'fusionmap':FusionMap,
    'jaffa':JAFFA,
    'mapsplice':MapSplice,
    'nfuse':nFuse,
    'soapfuse':SOAPfuse,
    'tophatfusion':TopHatFusion,
    'starfusion':STARFusion
}
