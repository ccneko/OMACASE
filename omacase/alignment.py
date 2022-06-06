#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20210213


"""
Handle OM alignment IO (currently read input only)
By Claire Chung
"""


import logging
from pathlib import Path
from collections import OrderedDict
from itertools import groupby
from operator import itemgetter

#import modin.pandas as pd
from intervaltree import Interval, IntervalTree

from omacase.opmap import Opmaps, BNX, CMAP, get_opmaps, segments_to_labels, pd

# ../m3e/SKBR3-DLE1/contigs/exp_refineFinal1_sv/EXP_REFINEFINAL1_full.xmap
logger = logging.getLogger(__name__)


class Cigar(object):
    """
    Cigar class to read alignment CIGAR strings;
    Modified from a short script in https://github.com/brentp/cigar 
    for sequence read CIGAR strings, which is under MIT license
    """
    q_consuming_ops = ("M", "I")
    r_consuming_ops = ("M", "D")

    def __init__(self, cigar_string):
        self.cigar = cigar_string

    def items(self):
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int(''.join(n)), ''.join(next(cig_iter)[1])

    def __str__(self):
        return self.cigar

    def __repr__(self):
        return "Cigar('%s')" % self

    def __len__(self):
        return sum(l for l, op,in self.items() if op in Cigar.q_consuming_ops)

    def r_len(self):
        return sum(l for l, op in self.items() if op in Cigar.r_consuming_ops)

    def reverse_cigar(self):
        return Cigar.string_from_elements(list(self.items())[::-1])
    
    def merge_like_ops(self):
        cigs = []
        for op, grps in groupby(self.items(), itemgetter(1)):
            cigs.append((sum(g[0] for g in grps), op))
        return Cigar(self.string_from_elements(cigs))

    @classmethod
    def string_from_elements(self, elements):
        return "".join("%i%s" % (l, op) for l, op in elements if l!=0)


class OMAlignment(object):
    """
    Optical Mapping alignment data I/O
    """
    def __init__(   self, align_id=None, r_id=None, q_id=None,
                    r_startpos=None, q_startpos=None,
                    r_endpos=None, q_endpos=None,
                    r_startseg=None, q_startseg=None,
                    r_endseg=None, q_endseg=None,
                    strand=None, cigar=None, confidence=None,
                    scheme='bng'):
        self.align_id    = None # only for XMAP
        self.r_id       = None
        self.q_id       = None
        self.r_startpos = -1
        self.q_startpos = -1
        self.r_endpos   = -1
        self.q_endpos   = -1
        self.r_startseg = -1
        self.q_startseg = -1
        self.r_endseg   = -1
        self.q_endseg   = -1
        self.strand     = None
        self.cigar      = None
        self.confidence = -1
        self.scheme     = None # {'bng', 'omtools'}

    @property
    def r_aligned_length(self):
        return abs(self.r_endpos - self.r_startpos)

    @property
    def q_aligned_length(self):
        return abs(self.q_endpos - self.q_startpos)

    @property
    def r_aligned_seg_num(self):
        return abs(self.r_endseg - self.q_endseg)

    @property
    def q_aligned_seg_num(self):
        return abs(self.q_endseg - self.q_endseg)

    @property
    def r_labels(self):
        return

    @property
    def q_labels(self):
        return

    @property
    def r_segments_aligned(self):
        return

    @property
    def q_segments_aligned(self):
        return

    @property
    def r_segments_unaligned(self):
        return

    @property
    def q_segments_unaligned(self):
        return

    @property
    def r_labels_aligned(self):
        return

    @property
    def q_labels_aligned(self):
        return

    @property
    def r_labels_unaligned(self):
        return

    @property
    def q_labels_unaligned(self):
        return


class OMAlignments(object):
    """
    OM alignment dataset; Collection of OMAlignment objects
    """
    def __init__(   self,
                    r_opmaps       =   None,
                    q_opmaps       =   None,
                    ivtree         =   None):
        self.r_opmaps       =   r_opmaps
        self.q_opmaps       =   q_opmaps
        self.ivtree         =   ivtree
        # sort by chrom, startpos, alignedlen, then qid, give alignid
        if self.ivtree is not None:
            self.get_ivtree()

    def to_oma(self):
        OMA().write_oma()

    def to_xmap(self):
        XMAP().write_oma()

    def stat(self):
        return self

    def filter_by_confidence(self):
        return self

    def filter_by_coverage(self):
        return self


class XMAP(OMAlignments):
    """
    XMAP Bionano map alignment data class;
    format specs from Bionano Document 30040, Rev. B
    """
    default_columns     =   [
                                "XmapEntryID", "QryContigID", "RefContigID",
                                "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                                "Orientation", "Confidence", "HitEnum",
                                "QryLen", "RefLen", "LabelChannel", "Alignment"
                            ]
    default_coltypes    =   [
                                int, int, int, float, float,
                                float, float, str, float, str,
                                float, float, int, str
                            ]
    def __init__(   self,
                    args        =   None,
                    filepath    =   None,
                    r_filepath  =   None,
                    q_filepath  =   None):
        self.filepath           =   Path(filepath)
        print(self.filepath)
        assert self.filepath.suffix.lower() == '.xmap'
        if r_filepath is not None:
            self.r_filepath     =   Path(r_filepath)
            self.r_opmaps       =   get_opmaps(self.r_opmaps)
        if q_filepath is not None:
            self.q_filepath     =   Path(q_filepath)
            self.q_opmaps       =   get_opmaps(self.q_opmaps)
        self.file_version       =   0.2
        self.filepath           =   None
        self.header_lines       =   None
        
        self.default_meta       =   OrderedDict()
        self.header             =   ''
        self.extra_comments     =   []
        self.strand             =   None
        self.rlabels           =   None
        self.qlabels           =   None
        self.read_oma()

    def load_ref(self):
        """Load CMAP reference"""
        if self.r_filepath is not None:
            self.r_opmaps   =   CMAP(self.r_filepath)
        elif self.r_opmaps is not None:
            assert isinstance(self.r_opmaps, Opmaps)
        else:
            logger.error('Invalid reference input. Please try again.')
        return self

    def load_query(self):
        """Load query maps"""
        if self.q_filepath is not None:
            self.q_opmaps   =   CMAP(self.q_filepath)
        elif self.q_opmaps is not None:
            assert isinstance(self.q_opmaps, Opmaps)
        else:
            logger.error('Invalid reference input. Please try again.')
        return self

    def get_ivtree(self):
        """Get alignment interval tree"""
        self.load_ref()
        self.load_query()
        self.ivtree = IntervalTree()
        return self

    def read_oma(self):
        f = open(self.filepath)
        for line in f:
            if line[0] == '#':
                self.header_lines.append(line)
            else:
                entry           =   line.split()
                align_id        =   int(entry[0])
                q_id            =   int(entry[1])
                r_id            =   int(entry[2])
                q_startpos      =   float(entry[3])
                q_endpos        =   float(entry[4])
                r_startpos      =   float(entry[5])
                r_endpos        =   float(entry[6])
                strand          =   entry[7]
                conf            =   float(entry[8])
                cigar           =   entry[9]
                q_len           =   entry[10]
                r_len           =   entry[11]
                labelchannel    =   entry[12]
                alignment       =   entry[13].strip('(').strip(')').split(')(')
                print(alignment[0])
                # execute this or cigar
                align           =   OMAlignment(
                                        align_id=align_id, r_id=r_id, q_id=q_id,
                                        r_startpos=r_startpos, q_startpos=q_startpos,
                                        r_startseg = alignment[0][0],
                                        q_startseg = alignment[0][1], strand=strand,
                                        cigar=cigar, confidence=conf, scheme='bng'
                                        )
                yield align
                # just need to read in cigar and first alignment
                # maybe assert with last alignment & total number

        f.close()
        return self

    def solve_cigar(self):
        return self

    def solve_alignment_string(self, alignment_string):
        return self

    def write_oma(self, filepath='output.xmap'):
        filepath = Path(filepath)
        assert str(self.filepath.suffix).lower() == '.xmap'
        if Path(self.filepath).is_file():
            raise FileExistsError
        f = open((self.filepath), 'w')
        rows = ''  # placeholder
        f.write(self.header+'\n')
        f.write(rows)
        f.close()


class OMA(OMAlignments):
    """
    OMA data class; child of OMAlignments class
    """
    def __init__(   self,
                    filepath    =   None,
                    r_filepath  =   None,
                    q_filepath  =   None):
        """Initializes an OMA object"""
        super().__init__()
        self.filepath           =   Path(filepath)
        self.r_filepath         =   Path(r_filepath)
        self.q_filepath         =   Path(q_filepath)
        assert str(self.filepath.suffix).lower() == '.oma'
        self.file_version       =   -1 #1.1
        self.extra_comments     =   None
        self.default_columns    =   [
                                        'QueryID', 'QuerySeg', 'QuerySegInfo', 'RefID',
                                        'Strand', 'Score', 'Confidence',
                                        'RefSegStart', 'RefSegStop', 'QuerySegStart',
                                        'RefStartCoord', 'RefStopCoord', 'Cigar'
                                    ]
        self.columns            = self.default_columns
        self.read_oma()

    def header(self):
        comments = self.extra_comments + '#' + '\t'.join()
        return

    def to_df(self, filepath=''):
        filepath = Path(filepath)
        assert str(filepath.suffix).lower() == '.oma'
        df = pd.read_csv(filepath, sep='\t', comment="#")
        df.columns = self.columns
        return df

    def read_oma(self, filepath=''):
        filepath = Path(filepath)
        assert str(filepath.suffix).lower() == '.oma'
        f = open(filepath)
        for line in f:
            if line[0] == '#':
                self.extra_comments += line
            else:
                break
        f.close()
        df = pd.read_csv(filepath, sep="\t", comment="#")
        for row in df.iterrows():
            yield OMAlignment()
        return

    def write_oma(self, filepath=None):
        filepath = Path(self.filepath)
        assert str(filepath.suffix).lower() == '.oma'
        if filepath.is_file():
            raise FileExistsError


def read_compnt_keys(filepath):
    keys = pd.read_csv(filepath, comment='#', sep='\t')
    keys = keys.set_index('CompntId', drop=True)
    return keys['CompntName'].to_dict()

def add_query_match_coords_to_oma_df(oma_df):
    oma_df['QueryStartCoord'] = oma_df.apply(lambda row: segments_to_labels(list(map(int, row['QuerySegInfo'].split(';'))))[row['QuerySegStart']], axis=1)
    oma_df['QueryStopCoord'] = oma_df.apply(lambda row: segments_to_labels(list(map(int, row['QuerySegInfo'].split(';'))))[row['QuerySegStop']], axis=1)
    return oma_df


def count_cigar(cigar):
    q_seg_cnt = 0
    r_seg_cnt = 0
    for i in range(0,len(cigar),2):
        if cigar[i+1] == 'M':
            q_seg_cnt += int(cigar[i])
            r_seg_cnt += int(cigar[i])
        elif cigar[i+1] == 'I':
            q_seg_cnt += int(cigar[i])
        elif cigar[i+1] == 'D':
            r_seg_cnt += int(cigar[i])
    return q_seg_cnt, r_seg_cnt