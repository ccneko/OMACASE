#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220219


"""
Structural Variations (SV) data I/O
"""


import logging
from pathlib import Path
from collections import defaultdict

import pandas as pd
from intervaltree import Interval, IntervalTree

from omacase.opmap import set_threads, get_opmaps, Opmaps
from omacase.utils import check_file_extension

logger = logging.getLogger(__name__)


def get_anno(args=None, filepath=None, from_string=False, string=None):
    """OMACASE main function call"""
    if args is not None:
        set_threads(args.threads)
        if args.input_file:
            logger.debug(args.input_file)
            filepath    =   args.input_file
    file_format, compressed =   check_file_extension(filepath, datatype = 'annotation')
    if file_format == 'bed':
        anno  =   BED(  filepath=filepath, from_string=from_string,
                        string=string, compressed=compressed)
    elif file_format == 'gtf':
        anno  =   GTF(  filepath=filepath, from_string=from_string,
                        string=string, compressed=compressed)
    else:
        raise ValueError('Please input filepath for supported annotation files.')
    logger.info(f'{len(anno)} annotations read')
    return anno


class Annotations(object):
    def __init__(self, ivtree=None):
        self.df             =   None
        self.ivtrees        =   defaultdict()
        self.colors         =   None

    def __len__(self):
        return len(self.df)


class BED(Annotations):
    def __init__(   self, filepath=None, threads=1,
                    header=None, columns_to_read=3,
                    from_string=False, string=False, compressed=False):
        super().__init__()
        self.filepath           =   Path(filepath)
        self.col_num            =   columns_to_read
        self.header             =   header
        self.read_anno() # set df & ivtrees
        self.__len__            =   len(self.df)
    
    def read_anno(self, data_colnames=None):
        self.df = pd.read_csv(self.filepath, sep='\t', comment='#', header=self.header, index_col=None)
        if data_colnames is None:
            self.df.columns     =   ['frag', 'start', 'end'] + list(self.df.columns[3:])
        else:
            self.df.columns     =   ['frag', 'start', 'end'] + data_colnames
        return self

    def get_ivtrees(self):
        for frag in self.df['frag'].unique():
            self.ivtrees[frag]  =   IntervalTree()
        if self.col_num == 3:
            for i, row in self.df.iterrows():
                self.ivtrees[row['frag']][row['start']:row['end']] = f'Feature {i}'
        elif self.col_num > 3:
            if self.col_num == 4 and str(self.df.columns[3]) == '3':
                self.df.columns     =   ['frag', 'start', 'end', 'data']
            for i, row in self.df.iterrows():
                self.ivtrees[row['frag']][row['start']:row['end']] = dict(row[3:])
        return self


class GTF(Annotations):
    """Initializeds an GTF object; child class of Annotations"""
    header  =   [   'seqname', 'source', 'feature', 'start', 'end',
                    'score', 'strand', 'frame', 'attribute' ]
    def __init__(   self, filepath=None, attributes=None, gencode=False,
                    from_string=False, string=False, compressed=False):
        super().__init__()
        self.filepath           =   Path(filepath)
        self.attributes         =   attributes
        self.read_anno(attributes=attributes, gencode=gencode) # set df & ivtrees

    def read_anno(self, attributes=None, gencode=False):
        self.df     =   pd.read_csv(    self.filepath, sep='\t', comment='#',
                                        header=None, index_col=None)
        self.df.columns = self.header
        if gencode  ==   True:
            self.attributes =   [   'gene_id', 'transcript_id', 'gene_type', 'gene_name', \
                                    'transcript_type', 'transcript_name', 'level', \
                                    'transcript_support_level', 'tag', 'havana_gene', \
                                    'havana_transcript' ]
        if self.attributes is not None:
            for attrib in self.attributes:
                self.df.loc[:, attrib]  =   self.df['attribute'].str.extract(\
                                                attrib + ' \"(\w+)\"', expand=False)
        return self
    
    def get_ivtrees(self):
        for i, row in self.df.iterrows():
            if row['seqname'] not in self.ivtrees:
                self.ivtrees[row['seqname']] = IntervalTree()
            self.ivtrees[row['seqname']][row['start']:row['end']] = dict(row[3:])
        return self.ivtrees


class SV(object):
    """Initailizes an SV object"""
    def __init__(self):
        self.entry_id       =   -1
        self.size_change    =   0
        self.type           =   None # insertion, deletion, inversion, translocation
        self.zygosity       =   None # homozygous, heterozygous, unknown
        self.score          =   -1
        self.scoring_scheme =   None # OMSV, Bionano


class SVs(Annotations):
    """
    Initializes an SVs object to store structural variation annotations;
    child class of Annotations
    """
    def __init__(self, filepath=None):
        super().__init__()
        self.filepath       =   Path(filepath)



class SMAP(SVs):
    """
    Bionano Document 30041, Rev. G
    https://bionanogenomics.com/wp-content/uploads/2017/03/30041-SMAP-File-Format-Specification-Sheet.pdf
    """
    supported_file_versions    =   ['0.8', '0.9']
    def __init__(self, filepath=None):
        super().__init__()
        self.filepath                   =   Path(filepath)
        self.file_version               =   -1
        
        self.columns    =   [   'SmapEntryID', 'QryContigID', 'RefcontigID1', 'RefcontigID2',
                                'QryStartPos', 'QryEndPos', 'RefStartPos', 'RefEndPos',
                                'Confidence', 'Type', 'XmapID1', 'XmapID2', 'LinkID',
                                'QryStartIdx', 'QryEndIdx', 'RefStartIdx', 'RefEndIdx',
                                'Zygosity', 'Genotype', 'GenotypeGroup','RawConfidence',
                                'RawConfidenceLeft', 'RawConfidenceRight', 'RawConfidenceCenter',
                                'SVsize', 'SVfreq', 'orientation' ]
        self.col_dtypes =   [   int, int, int, int, float, float, float, float, float, str,
                                int, int, int, int, int, int, int, str, int, int,
                                float, float, float, float, float, float, str ]
        self.read_omv()

    def read_omv(self):
        df          =   pd.read_csv(self.filepath, sep='\t', comment='#', header=None)
        df.columns  =   self.columns
        return self


class OSV(SVs):

    def __init__(self, filepath=None):
        super().__init__()
        self.filepath       =   Path(filepath)
        self.file_version   =   -1
        self.read_omv()

    def read_omv(self):
        return self


class VCF(Annotations):
    supported_file_versions =   ['4.0', '0.2']
    def __init__(self, filepath=None):
        super().__init__()
        self.filepath       =   Path(filepath)
        # 20190308
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        self.file_version   =   4.2
