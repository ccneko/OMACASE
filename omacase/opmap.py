#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220325


"""
Defines optical map data structures and IO in OMACASE
By Claire Chung
"""


import logging
import os, sys
import gc
import io
import re
import gzip
from pathlib import Path
from collections import OrderedDict, defaultdict
from multiprocessing import Manager, Pool, cpu_count
import math
import platform
import psutil

from intervaltree import Interval, IntervalTree
import numpy as np
from numba import jit
os.environ["MODIN_ENGINE"] = "ray"

import pandas as pd # avoid multiple ray init
from omacase.utils import log_func_run, assign_dtype, check_file_extension, timestamp, LOG_FORMAT


logging.basicConfig(format=LOG_FORMAT, stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)
logging.getLogger('numba').setLevel(logging.WARNING)


def set_threads(threads=cpu_count()):
    logger.info(f'Setting thread number to {threads}')
    os.environ["MODIN_CPUS"] = str(threads)
    import modin.pandas as pd
    import ray
    try:
        if platform.system() == 'Darwin':
            g=2
        else:
            g=4
        if not ray.is_initialized():
            ray.init(
                num_cpus=threads,
                object_store_memory=g*1000 * 1024 * 1024,
            )
    except OSError:
        raise

    return threads

def get_opmaps( args=None, filepath=None, string=None, filter_only=False,
                stat_only=False, labels_only=False, meta_only=False):
    """OMACASE main function call"""
    opmaps          =   None
    if args is not None:
        if hasattr(args, 'verbose') and args.verbose:
            logging.basicConfig(format=LOG_FORMAT, stream=sys.stderr, level=logging.DEBUG)
        if hasattr(args, 'input_file') and args.input_file:
            filepath    =   args.input_file
    if args is not None and hasattr(args, 'threads') and args.threads:
            threads = set_threads(args.threads)
    else:
        threads = set_threads()
    file_format, compressed =   check_file_extension(filepath)
    if file_format == 'bnx':
        opmaps  =   BNX(    
                            filepath=filepath,
                            string=string,
                            compressed=compressed,
                            stat_only=stat_only,
                            labels_only=labels_only,
                            meta_only=meta_only,
                            threads=threads
                        )
    elif file_format == 'cmap':
        opmaps  =   CMAP(   
                            filepath=filepath,
                            string=string,
                            compressed=compressed,
                            stat_only=stat_only,
                            labels_only=labels_only
                        )
    elif file_format == 'data':
        opmaps  =   OMDATA( 
                            filepath=filepath,
                            string=string,
                            compressed=compressed,
                    )
    else:
        ValueError('Please input filepath for supported OM dataset files.')
    if opmaps is None:
        logger.warning('No OM dataset was imported. Please import a dataset to proceed.')
    elif not isinstance(opmaps, Opmaps):
        raise ValueError(
            f'Input (type {type(opmaps)})' + \
            'is not Opmaps OM dataset object.')
    return opmaps


class Opmap(object):
    """
    Optical map data class object
    """
    def __init__(self, omid=0, segments=None, labels=None,
                 extra_info=None):
        """Initializes an Opmap object"""
        self.omid               =   omid
        self.segments           =   segments
        self.labels             =   labels
        self.extra_info         =   extra_info
        self.file_format        =   None
        self.label_density      =   None
        self.label_snr          =   None
        self.mol_snr            =   None
        self.label_intensity    =   None
        self.mol_intensity      =   None

        assert isinstance(segments, (list, type(None)))
        if segments is not None:
            assert isinstance(segments[0], (int, float))
            self.segnum     =   len(self.segments)
            self.labels     =   self.get_labels()

        if labels is not None:
            try:
                assert isinstance(labels[0], (int, float))
            except:
                logger.debug(labels)
                logger.critical(f'labels[0] {labels[0]} is not a number')
                raise

    def get_labels(self):
        """Get Opmap label label coordinates"""
        label_list          =   []
        pointer             =   0
        for label in self.segments:
            pointer         =   pointer + label + 1
            label_list.append(pointer)
        return label_list

    def get_segments(self, use_int=True):
        """Get Opmap segment lengths"""
        return labels_to_segments(self.labels, use_int=use_int, start=True)

    def __len__(self):
        """
        Length in bp. Call len(Opmap.segments) for segment number or 
        len(Opmap.labels) for label number.
        """
        #assert self.labels[-1] == sum(self.segments) + len(self.segments) - 1
        if self.labels is None:
            logger.warning(f'Molecule {self.omid} has zero length. Possibly an incomplete record')
            return 0
        else:
            return int(self.labels[-1])

    def __repr__(self):
        """Report an Opmap object as string in print()"""
        om_str  =   f'Map ID:\t{self.omid}\nMap length:\t{len(self)}\n'
        return om_str

    def __iter__(self):
        """Iterate segments in an Opmap object"""
        return iter(self.get_segments())

    def __add__(self, opmap_to_add):
        """Sum Opmap objects by joining"""
        return self.join(opmap_to_add)

    def __truediv__(self, divisor):
        """Truediv: Split Opmap per length"""
        return self.split_per_len(len(self) / divisor)

    def __floordiv__(self, divisor):
        """Floordiv: split Opmap by length"""
        return self.split_per_len(len(self) / divisor)

    def __getitem__(self, index):
        """Subset by start & end coordinates."""
        return self.subset_by_segments(index)

    @classmethod
    def from_seq(cls, seq, enzyme):
        """In silico digestion to create Opmap from sequence; placeholder"""
        opmap = seq.to_opmap() # placeholder
        return opmap

    def subset_by_segments(self, index):
        """Subset Opmap by segment index"""
        if self.segments is None:
            self.segments = self.get_segments()
        logger.debug(f'Extracting segments: {index}')
        new_segments: list  =   self.segments[index]
        if isinstance(index, int):
            return self.segments[index]
        if new_segments[0] == 19:
            new_segments.pop(0)
        if new_segments[-1] == 20:
            new_segments    =   new_segments[:-1]
        start_index         =   index.start
        end_index           =   index.stop
        if start_index is None:
            start_index     =   0
        elif start_index < 0:
            start_index     =   len(self.segments) + index.start
        if end_index is None:
            end_index       =   len(self.segments) - 1
        elif end_index < 0:
            end_index       =   len(self.segments) + index.stop

        new_id              =   str(self.omid) + '_seg' + str(start_index) + '_seg' + str(end_index)
        new_opmap           =   Opmap(new_id, new_segments)
        return new_opmap

    def extract_by_segment(self, start_seg, end_seg=-1):
        """Extract subsection of an Opmap by segment index"""
        logger.info('0-based from starting segment 19.')
        logger.info('Input label number to be split at.')
        if self.segments is None:
            self.segments = self.get_segments()
        assert end_seg >= 0 and end_seg < len(self.segments)
        if start_seg == 0:
            start_seg       =   1
        if end_seg == -1:
            end_seg         =   len(self.segments) - 1
        new_id              =   str(self.omid) + '_seg' + (str(start_seg)) + '_' + (str(end_seg))
        new_segments        =   self.segments[start_seg:end_seg + 1]
        assert start_seg != 0
        new_segments        =   [19] + new_segments
        if end_seg != len(self.segments) - 1:
            new_segments    =   new_segments + [20]
        new_opmap           =   Opmap(new_id, new_segments)
        return new_opmap

    def join(self, to_join_opmap, join_start_seg=1, join_end_seg=0, start_seg=1,
             end_seg=0, new_id='', how='right'):
        """Join 2 Opmaps. Subscriptable"""
        logger.info('Join 2 Opmaps. Under construction.')
        assert join_start_seg > 0
        assert join_end_seg >= 0 and end_seg < len(to_join_opmap.segments)
        assert start_seg > 0
        assert end_seg >= 0 and end_seg < len(self.segments)
        assert how == 'right' or how == 'left'
        if new_id == '':
            new_id           =   str(self.omid) + '_' + to_join_opmap.omid
        if end_seg == 0:
            end_seg         =   len(self.segments) - 2
        if join_end_seg == 0:
            join_end_seg    =   len(to_join_opmap.segments) - 2
        segfromself         =   self.segments[start_seg:end_seg + 1]
        segtoadd            =   to_join_opmap.segments[join_start_seg:join_end_seg + 1]
        logger.debug(to_join_opmap.segments[start_seg:end_seg + 1])
        logger.debug(segtoadd[join_start_seg:join_end_seg + 1])
        new_segments = []
        if how == 'right':
            new_segments    =   segfromself + segtoadd
        else:
            new_segments    =   segtoadd + segfromself
        new_segments        =   [19] + new_segments + [20]
        new_opmap           =   Opmap(new_id, new_segments)
        return new_opmap

    def flip(self):
        """Flip Opmap direction"""
        if int(float(self.segments[0])) == 19 and \
            int(float(self.segments[-1])) == 20:
            new_segments = [19]+self.segments[1:-1][::-1]+[20]
        else:
            new_segments        =   self.segments[::-1]
        new_opmap           =   Opmap(self.omid + 1e9, new_segments)
        return new_opmap

    def scale(self, scale_factor=1):
        """Scale Opmap"""
        new_segments        =   list(map(lambda x: int(x * scale_factor), self.segments))
        new_segments[0]     =   19
        new_segments[-1]    =   20
        return Opmap(self.omid, new_segments)

    def ivtree(self):
        """Create an intervaltree object from Opmap"""
        tree                =   IntervalTree(Interval(  self.labels[i], self.labels[i+1])
                                                        for i in range(len(self.labels)-1))
        return tree

    def extract_by_len(self, start_pos=0, end_pos=0):
        """Extract subsection of an Opmap by coordinates"""
        assert end_pos > start_pos
        new_id                  =   str(self.omid) + '_' + str(start_pos) + '_' + str(end_pos)
        t                       =   self.ivtree()
        ivs                     =   sorted(list(t[start_pos:end_pos]))
        if len(ivs) == 0:
            new_segments        =   [end_pos - start_pos]
            new_opmap           =   Opmap(new_id, new_segments)
            return new_opmap
        startiv                 =   ivs[0]
        endiv                   =   ivs[-1]
        logger.debug(startiv.begin, endiv.begin)
        startidx                =   sorted(list(t)).index(startiv)
        endidx                  =   sorted(list(t)).index(endiv) + 2
        if endidx == startidx:
            endidx              +=  1
        new_segments            = self.segments[startidx:endidx]
        new_segments[0]         -=  start_pos - sum(self.segments[:startidx])
        if new_segments[-1] > abs(end_pos - sum(self.segments[:endidx])):
            new_segments[-1]    -=  abs(end_pos - sum(self.segments[:endidx]))
        if sum(new_segments) > end_pos - start_pos:
            pass
        new_opmap               =   Opmap(new_id, new_segments)
        check                   =   sum(new_opmap.segments) - len(new_opmap) \
                                    + len(new_opmap.segments) - 1
        assert check == 0
        return new_opmap

    def split_per_len(self, divisor):
        """Split an Opmap into multiple Opmap objects per provided length); Under construction"""
        new_opmaps = []
        for i in range(divisor, len(self), divisor):
            logger.debug([i - divisor, i - 1])
            new_opmaps.append(self.extract_by_len(i - divisor, i - 1))
        return new_opmaps

    def head(self, n=10):
        return self[:n]

    def tail(self, n=10):
        return self[-n:]


class Opmaps(object):
    """
    OM dataset data class, collection of Opmap objects
    """
    def __init__(self):
        self.omdict             =   defaultdict(Opmap)
        self.filepath           =   None
        self.filename           =   None
        self.mapnum             =   None
        self.stat_only          =   False
        self.labels_only        =   False
        self.df                 =   None
        self.lengths            =   []
        self.label_density      =   0
        self.n50                =   None
        self.n90                =   None

    def get_lengths(self):
        """Get sum of Opmap lengths in a CMAP dataset object"""
        if len(self.lengths) > 0:
            return self.lengths
        elif self.df is None:
            return np.array([len(self.omdict[omid]) for omid in self.omdict])
        elif self.df is not None:
            return self.df['length'].to_numpy()
    
    @property
    def total_length(self):
        """Get number of Opmap objects in an Opmaps dataset object"""
        return sum(self.lengths)

    def __len__(self):
        """Get number of Opmap objects in an Opmaps dataset object"""
        if self.df is not None:
            return len(self.df)
        else:
            return len(self.omdict)
        
    def __getitem__(self, index):
        """Get Opmap objects or df rows from an Opmaps dataset object by index"""
        if self.df is not None:
            return self.df.loc[:, index]
        else:
            return self.omdict[index]

    @classmethod
    def from_fasta(cls, filepath, enzyme=None):
        return FASTA.to_opmaps(filepath, enzyme)

    def to_bnx(self, filepath, file_version):
        """Output Opmaps object in BNX format"""
        BNX(omdict=self.omdict).write_om(filepath=filepath, file_version=file_version)

    def to_cmap(self, filepath, file_version):
        """Output Opmaps object in CMAP format"""
        #CMAP(omdict=self.omdict).write_om(filepath=filepath, file_version=file_version)
        raise NotImplementedError

    def to_omdata(self, filepath, file_version):
        """Output Opmaps object in OMDATA format"""
        #OMDATA(omdict=self.omdict).write_om(filepath=filepath, file_version=file_version)
        raise NotImplementedError

    def rm_opmaps_by_ids(self, to_rm):
        """Remove an entry from OM dataset by ID"""
        [self.omdict.pop(omid) for omid in to_rm]

    def summary(self):
        """Report OM dataset statistics"""
        logger.info(f'Number of maps:\t{len(self)}')
        logger.info(f'Total length:\t{self.total_length}')
        logger.info(f'N50:\t{self.get_n50()}')
        logger.info(f'label density:{self.label_density}')
        return

    def n_metrics(self, n=50):
        """Calculate and return N-metrics, by default N50, for an input omdict"""
        return calc_n_metrics(self.lengths, n)

    @property
    def keys(self):
        """Returns keys of Opmap objects in omdict in subscriptable list"""
        return sorted(list(self.omdict.keys()))

    def get_n50(self):
        """Returns N50 and L50 values"""
        if self.n50 is not None:
            return self.n50
        else:
            self.n50 = self.n_metrics(50)
            return self.n50

    def get_n90(self):
        """Returns N90 and L90 values"""
        if self.n90 is not None:
            return self.n90
        else:
            self.n90 = self.n_metrics(90)
            return self.n90

    def get_label_density(self):
        """Returns average number of labels per 100 kbp"""
        if self.label_density is not None:
            return self.label_density
        else:
            if not self.stat_only:
                try:
                    label_num  =   sum([len(self.omdict[x].labels) for x in self.omdict])
                except TypeError:
                    logger.critical('TypeError exception occurred: tracing the erratic entry')
                    for x in self.omdict:
                        try:
                            len(self.omdict[x].labels)
                        except TypeError:
                            logger.exception(f'TypeError due to map ID {x}')
                    raise TypeError('Please check error log for details')
                return label_num / self.total_length * 100000

    def filter_by_label_density(self, min_val=-1, max_val=-1):
        """Filter members in an Opmaps dataset object by label density"""
        new_opmaps      =   Opmaps()
        if min_val == -1 and max_val == -1:
            logger.warning('No filtering criteria input. No filtering is done.')
            return self
        for opmap in self.omdict:
            if  (opmap.label_density >  min_val or min_val == -1) and \
                (opmap.label_density <= max_val or max_val == -1):
                new_opmaps.omdict[opmap.omid] = opmap.segments
        return new_opmaps

    def filter_by_length(self, min_val=-1, max_val=-1):
        """Filter members in an Opmaps dataset object by length"""
        new_opmaps      =   Opmaps()
        if min_val == -1 and max_val == -1:
            logger.warning('No filtering criteria input. No filtering is done.')
            return self
        for opmap in self.omdict:
            if  (len(opmap) >  min_val or min_val == -1) and \
                (len(opmap) <= max_val or max_val == -1):
                new_opmaps.omdict[opmap.omid] = opmap.segments
        return new_opmaps

    def get_segments(self, save_segments=False):
        """Get segment lengths per map from label positions"""
        if self.df is not None:
            labels = self.df['labels']._to_pandas()
            segments = labels.apply(labels_to_segments)
            if not save_segments:
                return segments
            else:
                self.df['segments'] = segments
                del segments


class OMDATA(Opmaps):
    """
    OMDATA data class; child of Opmap class
    """
    header  =   '#Fragment ID\tSize\tTotalSegments\tSegmentDetail'
    def __init__(   self, filepath=None, filename=None,
                    string=None, compressed=False,
                    stat_only=False, read_meta=True, read_map_num=-1):
        """Initializeds an OMDATA object"""
        super().__init__()
        self.read_meta          =   read_meta
        self.read_map_num       =   int(read_map_num)
        self.from_string        =   bool(string)
        if not self.from_string:
            self.filepath       =   filepath
            self.filename       =   Path(self.filepath).name
            logger.info(f'Start OMDATA file reading: {self.filepath}')
        else:
            self.filepath       =   None
            self.filename       =   Path(filename)
        self.extra_comments     =   []
        self.header             +=  '\n'.join(self.extra_comments)
        self.default_columns    =   ['frag_id', 'size', 'segments', 'segment_detail']
        self.read_om()
        logger.info(f'Complete OMDATA file reading.')

    @log_func_run
    def read_om(self):
        """Read OMDATA file into an OMDATA object"""
        assert Path(self.filename).suffix.lower() == '.data', \
        f'{Path(self.filename).suffix.lower()} format is not supported.'
        self.df                      =   pd.read_csv(self.filepath, sep='\t', comment='#', header=None)
        self.df.columns              =   self.default_columns
        self.df.set_index('frag_id')
        self.df['segment_detail']    =   self.df['segment_detail'].str.split(';')
        self.df['segment_detail']    =   self.df['segment_detail'].apply(lambda x: list(map(int, x)))
        self.omdict = {}
        for row in self.df.iterrows():
            self.omdict[row[0]] =   Opmap(omid=row[0], segments=row[1]['segment_detail'])
        return self

    @log_func_run
    def write_om(self, omids=None, filepath=sys.stdout, file_version='0.2'):
        """Write OMDATA file output"""
        if filepath != sys.stdout: # to simplify with decorator
            assert Path(filepath).suffix.lower in ['.data']
        elif Path(filepath).is_file():
            raise FileExistsError
        f = open(self.filepath, 'w')
        f.write(self.header+'\n')
        f.write('\n'.join([str(self.omdict[x]) for x in self.omdict]))
        f.close()


class BNX(Opmaps):
    """
    BNX data class; child of Opmap class; format specs from Bionano Document 30038, Rev. F
    """
    supported_file_versions =   ['1.2', '1.3']
    row0_columns_common     =    [
                                    'LabelChannel', 'MoleculeId', 'Length',
                                    'AvgIntensity', 'SNR', 'NumberofLabels',
                                    'OriginalMoleculeId', 'ScanNumber',
                                    'ScanDirection', 'ChipId', 'Flowcell', 'RunId'
                                ]
    row0_columns_by_ver     =   {
                                    '1.2': ['GlobalScanNumber'],
                                    '1.3': ['Column', 'StartFOV', 'StartX',
                                            'StartY', 'EndFOV', 'EndX', 'EndY',
                                            'GlobalScanNumber']
                                }
    row0_types_common       =   [
                                    int, int, float, float, float,
                                    int, int, int, int, str,
                                    int, int
                                ]
    row0_types_by_ver       =   {
                                    '1.2': [int],
                                    '1.3': [int] * 8
                                }
    col_dtype_desc          =   []
    # Data rows
    rown_columns           =   ['LabelChannel', 'LabelPositions[n]']
    rown_types             =   [str, float]
    row0_types             =   None
    # QX11/QX21 for label SNR, QX12/QX22 for label intensity
    rowqxnn_columns        =   ['QualityScoreID', 'QualityScores[n]']
    rowqxnn_types          =   [str, float]
    extra_comments         =   []
    
    def __init__(   self, filepath=None, filename=None,
                    string=None, compressed=False, omdict=None, 
                    stat_only=False, labels_only=False, 
                    meta_only=False, filter_only=False,
                    read_meta=True, read_map_num=-1, read_extra_cols=False,
                    threads=cpu_count(),
                    ):
        """Initializes a BNX object"""
        super().__init__()
        self.threads                    =   threads
        self.file_format                =   'bnx'
        if omdict is not None:
            self.omdict                 =   omdict
        self.df                         =   None
        self.meta_only                  =   meta_only
        self.filter_only                =   filter_only
        self.labels_only                =   labels_only
        self.stat_only                  =   stat_only
        self.approx_stat                =   True
        self.read_meta                  =   read_meta
        self.read_map_num               =   int(read_map_num)
        self.from_string                =   bool(string)
        if not self.from_string:
            self.filepath               =   Path(filepath)
            self.filename               =   Path(self.filepath).name
        else:
            self.filepath               =   None
            self.filename               =   Path(filename)
        self.compressed                 =   compressed
        self.read_extra_cols            =   read_extra_cols
        self.file_version               =   '0.0'
        self.label_channels             =   1
        self.motifs                     =   []
        self.colors                     =   []
        self.software_version           =   None # optional
        self.bases_per_pixel            =   None # optional
        self.default_meta               =   OrderedDict()
        self.chips                      =   []
        self.qx11_avg                   =   -1
        self.qx12_avg                   =   -1
        self.qx21_avg                   =   -1
        self.qx22_avg                   =   -1
        self.mol_snr_avg                =   -1
        self.mol_intensity_avg          =   -1
        self.qx11_sd                    =   -1
        self.qx12_sd                    =   -1
        self.qx21_sd                    =   -1
        self.qx22_sd                    =   -1
        self.mol_snr_sd                 =   -1
        self.mol_intensity_sd           =   -1
        
        if self.from_string:
            self.read_om(
                string=string,
            )
        else:
            self.read_om()
        if not meta_only:
            if not self.stat_only:
                self.mapnum             =   len(self.omdict)
            else:
                self.mapnum             =   len(self.df)
            if self.mapnum is not None:
                if len(self.omdict) != self.mapnum  and self.mapnum != len(self.df):
                    logger.warning('Total number of molecules read does not match with description.')
            logger.debug(f'lengths nan num: {np.count_nonzero(np.isnan(self.lengths))}')
            self.get_label_density()

    @log_func_run
    def extract_meta(self, lines):
        """Extract metadata from # comments in a BNX file"""
        logger.info(f'Start reading BNX file metadata')
        file_ver_line_pattern               =   re.compile(r'[Ff]ile [Vv]ersion:[\t\s]+(\d+\.\d+)')
        label_channels_line_pattern         =   re.compile(r'[Ll]abel [Cc]hannels:[\t\s]+(\d+)')
        map_num_line_pattern                =   re.compile(r'[Nn]umber of [Cc]onsensus [Mm]aps:[\t\s]+(\d+)')
        logger.debug(f'Header line number: {len(lines)}')
        for i, line in enumerate(lines):
            if line.split('\t', 1)[0] == '#0h':
                self.row0_columns = line.split('\t', 1)[1].split()
                if len(self.row0_columns) >= 20:
                    self.file_format = '1.3'
                else:
                    self.file_format = '1.2'
                continue
            if  self.file_version   ==  '0.0':
                match_file_ver              =    re.search(file_ver_line_pattern, line[1:].strip())
                if  match_file_ver  is not None:
                    self.file_version       =   match_file_ver[1]
                    continue
            match_channels              =   re.search(label_channels_line_pattern, line[1:].strip())
            if match_channels   is not None:
                self.label_channels     =   int(match_channels[1])
                continue
            match_recogn = None
            if len(self.motifs) < self.label_channels:
                match_recogn                =   re.search('[Nn]ickase [Rr]ecognition [Ss]ite \d' \
                                                + ':[\t\s]+(\S+)\;?(\S+)?', line )
                if match_recogn is not None:
                    if match_recogn[1]  is not None:
                        self.motifs.append(match_recogn[1])
                    if match_recogn[2]  is not None:
                        self.colors.append(match_recogn[2])
                    continue
            match_map_num                   =   re.search(map_num_line_pattern, line)
            if  match_map_num       is not  None:
                self.mapnum                 =   int(match_map_num[1])
                continue
            if  match_file_ver      is None and match_channels  is None and \
                match_recogn        is None and match_map_num   is None:
                self.extra_comments.append(line.strip())
            if  line[:2] == '#h':
                self.columns                =   list(map(lambda x: x.strip(), line[3:].split('\t')))
                if len(self.columns) >= len(self.default_columns):
                    logger.warning('This BNX file contains non-standard columns.')
                    if not self.read_extra_cols:
                        logger.warning('Dropping info from extra columns.')
            if line[:2] == '#f':
                self.col_dtype_desc         =   line[2:].split('\t')
        if self.file_version not in self.supported_file_versions:
            logger.warning(   f'BNX file version {self.file_version} not supported.' +\
                                f'Supported versions are {self.supported_file_versions}')
        self.default_meta['file_version']   =   '# BNX File Version:\t' + self.file_version
        self.row0_columns                   =   self.row0_columns_common + \
                                                self.row0_columns_by_ver[self.file_version]
        logger.info(f'File format: {self.file_format.upper()} v{self.file_version}')
        self.row0_types                     =   self.row0_types_common \
                                                + self.row0_types_by_ver[self.file_version]

        self.default_meta['label_channels'] =   str(self.label_channels)
        self.default_meta['mapnum']         =   '# Number of Molecules:\t' + str(self.mapnum)
        for i in range(1,self.label_channels+1):
            if len(self.motifs) < i:
                self.motifs.append('unknown')
            if len(self.colors) < i:
                self.colors.append('unknown')
            self.default_meta['recogn_motif']   =   '# Nickase Recognition Site ' + \
                                                    str(i) + ':\t' + \
                                                    str(self.motifs[i-1]) + ';' + \
                                                    str(self.colors[i-1])
        self.default_meta['0h'] =   '#0h\t' + '\t'.join(self.row0_columns)
        self.default_meta['0f'] =   '#0f\t' + '\t'.join(self.col_dtype_desc)
        self.header             =   '\n'.join(self.extra_comments +
                                                    [self.default_meta[x] for x in self.default_meta])
        self.default_columns    =   self.row0_columns_common + \
                                    self.row0_columns_by_ver[self.file_version]
        return self

    def read_om(
            self, string=None,
        ):
        logger.info(f'Start BNX file reading: {self.filename}')
        if self.from_string:
            if self.compressed:
                string          =   gzip.decompress(string).decode('utf-8')
            f                   =   io.StringIO(string)
        else:
            if self.compressed:
                gz              =   gzip.open(self.filepath, 'r')
                f               =   io.BufferedReader(gz)
            else:
                f               =   open(self.filepath, 'r')
            file_bytes          =   os.stat(self.filepath).st_size
            bytes_per_chunk     =   math.ceil(file_bytes/self.threads)
            bytes_per_chunk     +=  bytes_per_chunk % 4
        header_lines = []
        lines = []
        while True:
            line = f.readline()
            if self.compressed:
                line = line.decode('utf-8')
            if line[0] == '#':
                header_lines.append(line[1:].strip())
            else:
                lines.append(line)
                break
        self.extract_meta(header_lines)
        if self.meta_only:
            return self
        logger.info('Collecting entry lines from file')
        lines                   +=  f.readlines()
        if isinstance(lines[1], bytes):
            lines = [lines[0]] + [line.decode('utf-8') for line in lines[1:]]
        lines                   = list(filter(lambda x: x[0] != '#', lines))
        lines                   = [str(math.floor(i/4)) + '\t' + line for i, line in enumerate(lines)]
        logger.debug(f'{lines[0][:10]}...')
        logger.info('Opening parallel thread pool for to read data rows')
        logger.info(f'Using {self.threads} threads')
        m = Manager()
        q = m.Queue()
        min_lines_to_read = 65536
        if 'cttaag' in ' '.join(self.motifs).lower() or 'dle' in ' '.join(self.motifs).lower():
            min_lines_to_read = 32768
        line_num_to_read        =   min(min_lines_to_read, int(len(lines)/self.threads))
        line_num_to_read        +=  4 - (line_num_to_read % 4)
        logger.debug(f'line_num_to_read: {line_num_to_read}')
        job_num                 =   math.ceil(len(lines)/line_num_to_read)
        logger.info(f'Total job number: {job_num}')
        input_iter  =   lambda i: ( 
                            i, lines[i*line_num_to_read:(i+1)*line_num_to_read],
                            self.row0_types, self.row0_columns, 
                            self.stat_only, q,
                        )
        pool        =   Pool(self.threads)     
        jobs        =   pool.starmap(   read_bnx_chunk_wrapper,
                                        [input_iter(i) for i in range(job_num)])
        assert len(jobs) == job_num
        
        if self.filter_only:
            jobs    =   pool.starmap(   filter_bnx_chunk_wrapper,
                                        [input_iter(i) for i in range(job_num)])
            [q.get() for job in jobs]
            return self

        tmp_bnxstat_lists = []
        for job in jobs:
            tmp_bnxstat_lists.append(q.get())
        pool.close()
        logger.info('All parallel jobs in queue complete. Pool closed.')
        logger.info('Converting to dataframe.')
        self.df = pd.DataFrame([item for sublist in tmp_bnxstat_lists for item in sublist]).dropna()
        del tmp_bnxstat_lists
        gc.collect()
        if self.stat_only:
            self.df.columns = [
                    'omid', 'length', 'label_num', 'mol_snr', 'mol_intensity',
                    'flowcell', 'scan', 'qx11_mol_mean', 'qx12_mol_mean'
            ]
        else:
            self.df.columns = [
                    'omid', 'length', 'label_num', 'mol_snr', 'mol_intensity',
                    'flowcell', 'scan', 'qx11_mol_mean', 'qx12_mol_mean', 'labels'
            ]
        logger.info('Complete converting to dataframe.')
        self.set_quality_stats(self.df)
        invalid_ids =   list(filter(lambda x: isinstance(x, str), self.omdict.keys()))
        valid_ids   =   list(filter(lambda x: isinstance(x, int), self.omdict.keys()))
        if len(invalid_ids) > 0:
            logger.warning(f'Number of invalid IDs: {len(invalid_ids)}')
            logger.warning(f'Invalid IDs: {invalid_ids}')

        if not self.stat_only and self.df is None:
            for x in invalid_ids:
                self.omdict[max(valid_ids) + 1] = self.omdict[x]
                del self.omdict[x]
                logger.warning(f'Invalid ID {x} changed to {max(valid_ids)+1}')
            self.omdict =   {k: self.omdict[k] for k in sorted(self.omdict)}

        logger.info(f'Complete BNX file reading: {self.filename}')
        self.lengths = self.get_lengths()
        return self

    def write_om(
            self, omids=None, 
            output_prefix=f'{timestamp()}-output', stdout=False,
            compressed=True, file_version='1.3'
        ):
        """Write BNX file output; on-going work"""
        if not stdout: # to simplify with decorator
            filepath = f'output_prefix' + 'compressed'*compressed
            if Path(filepath).is_file():
                raise FileExistsError

        f = open(filepath, 'w')
        bnx_rows =  ''  # placeholder, omdict
        header  =   '\n'.join(self.extra_comments +
                    [self.default_meta[x] for x in self.default_meta])
        f.write(header+'\n')
        f.write(bnx_rows)
        f.close()
        return self
    

    def get_label_density(self):
        """Get overall label density in a BNX dataset object"""
        if self.stat_only:
            self.label_density = sum(self.df['label_num'].to_numpy()) / sum(self.lengths) * 100000
        elif self.labels_only:
            self.label_density =  sum(self.df.loc[:,'labels'].apply(len))
        else:
            label_num  =   sum([len(self.omdict[x].labels) for x in self.omdict])
            self.label_density = label_num / self.total_length * 100000
        return self
    
    def filter_by_label_snr(self, min=-1, max=-1):
        """Filter BNX by SNR (default: min=-1, max=-1, i.e. not filtered)"""
        if min == -1 and max == -1:
            return self
        for omid in self.omdict.keys():
            if self.omdict[omid].label_snr < min or self.omdict[omid].label_snr > max:
                del self.omdict[omid]
        return self
    
    def filter_by_mol_snr(self, min=-1, max=-1):
        """Filter BNX by SNR (default: min=-1, max=-1, i.e. not filtered)"""
        if min == -1 and max == -1:
            return self
        for omid in self.omdict.keys():
            if self.omdict[omid].mol_snr < min or self.omdict[omid].mol_snr > max:
                del self.omdict[omid]
        return self

    def set_quality_stats(self, df):
        """Calculate QX11 (Channel 1 Label SNR) statistics"""
        if self.qx11_avg == -1 or self.qx11_sd == -1:
            qx11_all = df['qx11_mol_mean'].to_numpy()
            if self.qx11_avg == -1:
                self.qx11_avg = numba_mean(qx11_all)
            if self.qx11_sd == -1:
                self.qx11_sd = numba_std(qx11_all)

        """Calculate QX12 (Channel 1 Label intensity) statistics"""
        if self.qx12_avg == -1 or self.qx12_sd == -1:
            qx12_all = df['qx12_mol_mean'].to_numpy()
            if self.qx12_avg == -1:
                self.qx12_avg = numba_mean(qx12_all)
            if self.qx12_sd == -1:
                self.qx12_sd = numba_std(qx12_all)

        """Calculate QX21 (Channel 2 Label SNR) statistics"""
        if self.label_channels == 2 and (self.qx21_avg == -1 or self.qx21_sd == -1):
            qx21_all = df['qx21_mol_mean'].to_numpy()
            if self.qx21_avg == -1:
                self.qx21_avg = numba_mean(qx21_all)
            if self.qx21_sd == -1:
                self.qx21_sd = numba_std(qx21_all)

        """Calculate QX22 (Channel 2 Label intensity) statistics"""
        if self.label_channels == 2 and (self.qx22_avg == -1 or self.qx22_sd == -1):
            qx22_all = df['qx22_mol_mean'].to_numpy()
            if self.qx22_avg == -1:
                self.qx22_avg = numba_mean(qx22_all)
            if self.qx22_sd == -1:
                self.qx22_sd = numba_std(qx22_all)

        """Calculate molecule backbone SNR statistics"""
        if self.label_channels == 2 and self.mol_snr_avg == -1 or self.mol_snr_sd == -1:
            mol_snr_all = df['mol_snr'].to_numpy()
            if self.mol_snr_avg == -1:
                self.mol_snr_avg = numba_mean(mol_snr_all)
            if self.mol_snr_sd == -1:
                self.mol_snr_sd = numba_std(mol_snr_all)

        """Calculate molecule backbone intensity statistics"""
        if self.label_channels == 2 and self.mol_intensity_avg == -1 or self.mol_intensity_sd == -1:
            mol_intensity_all = df['mol_intensity'].to_numpy()
            if self.mol_intensity_avg == -1:
                self.mol_intensity_avg = numba_mean(mol_intensity_all)
            if self.mol_intensity_sd == -1:
                self.mol_intensity_sd = numba_std(mol_intensity_all)


def read_bnx_chunk_wrapper( job_id, lines_to_read,
                            row0_types, row0_columns,
                            stat_only, q):
    auto_garbage_collect()
    return q.put(read_bnx_chunk(
                    job_id, lines_to_read,
                    row0_types, row0_columns, stat_only))

def read_bnx_chunk( job_id, lines_to_read,
                    row0_types, row0_columns, stat_only):
    """Read BNX file into BNX object for statistics only"""
    tmp_bnxstat_list    =   []
    for i, line in enumerate(lines_to_read):
        if i!=0 and not i%10000:
            logger.debug(f'{line[:10]}...')
            logger.info(f'{i} molecules read in job {job_id}')
        split_line = line.split()
        omid = int(split_line[0])
        if split_line[1] == '0':
            tmp_bnxstat_list.append([])
            row0_values =   ['1'] + split_line[2:]
            row0_values =   [assign_dtype(x, (row0_types)[i]) \
                            for i, x in enumerate(row0_values)]
            mol_header  =   dict(zip(row0_columns, row0_values))
            mol_header['MoleculeId'] = omid
            chip = mol_header['ChipId'].strip('chips,')
            tmp_bnxstat_list[-1] += [
                mol_header['MoleculeId'],
                mol_header['Length'],
                mol_header['NumberofLabels'],
                mol_header['SNR'],
                mol_header['AvgIntensity'],
                mol_header['Flowcell'],
                mol_header['ScanNumber'],
                -1,
                -1
            ]
            if not stat_only:
                tmp_bnxstat_list.append([])
        elif split_line[1] == '1' and not stat_only:
            try:
                labels = list(map(float, split_line[2:]))
                tmp_bnxstat_list[-1][9] = labels
            except IndexError:
                tmp_bnxstat_list.append([omid]+[-1]*8)
                tmp_bnxstat_list[-1] += labels
        elif split_line[1][0] == 'Q':
            arr = np.array(list(map(float, split_line[2:])))
            if split_line[1] == 'QX11':
                if len(split_line) == 2:
                    qx11_mol_mean = 0
                else:
                    qx11_mol_mean = numba_mean(arr)
                try:
                    tmp_bnxstat_list[-1][7] = qx11_mol_mean
                except IndexError:
                    tmp_bnxstat_list.append([omid]+[-1]*6)
                    tmp_bnxstat_list[-1] += [qx11_mol_mean, -1]
            if split_line[1] == 'QX12':
                if len(split_line) == 2:
                    qx12_mol_mean = 0
                else:
                    qx12_mol_mean = numba_mean(arr)
                try:
                    tmp_bnxstat_list[-1][8] = qx12_mol_mean
                except IndexError:
                    tmp_bnxstat_list.append([omid]+[-1]*7)
                    tmp_bnxstat_list[-1] += [qx11_mol_mean]
        elif split_line[1][0] not in ['#', '0', '1', 'Q']:
                logger.warning(f'Invalid line detected: {line[:20]}...')
                logger.warning(f'A line in BNX file should start with #, 0, 1 or QX, but not {split_line[1]}')
    return tmp_bnxstat_list


def mol_snr_per_scan(opmaps):
    results = defaultdict(list)
    for key in list(opmaps.omdict.keys()):
        info = opmaps.omdict[key].extra_info
        results['_'.join(list(map(str,[info['ChipId'],info['Flowcell'],info['RunId']])))]\
        .append(opmaps.omdict[key].extra_info['SNR'])
    for k, v in results.items():
        results[k] = sum(v)/len(v)
    return results


class CMAP(Opmaps):
    """
    CMAP data class; child of Opmap class; format specs from Bionano Document 30039, Rev. G
    """
    supported_file_versions = ['0.1', '0.2']
    default_columns         =   [
                                        'CMapId', 'ContigLength', 'NumSites',
                                        'SiteID', 'LabelChannel', 'Position',
                                        'StdDev', 'Coverage', 'Occurrence',
                                        'ChimQuality', 'SegDupL', 'SegDupR',
                                        'FragileL', 'FragileR', 'OutlierFrac',
                                        'ChimNorm', 'Mask',
                                ]
    col_dtypes         =   [int, float, int, int, int] + \
                                    [float] * 11 + [str] # last str is hex
    col_dtype_desc     =   ['int', 'float'] + ['int']*3+ ['float']*11 + ['Hex']
    def __init__(   self, filepath=None, filename=None,
                    string=None,
                    compressed=False, omdict=None,
                    stat_only=False, labels_only=True,
                    read_meta=True, read_map_num=-1, read_extra_cols=False):
        """Initializes a CMAP object"""
        super().__init__()
        self.file_format        =   'cmap'
        self.from_string        =   bool(string)
        if omdict is not None:
            self.omdict         =   omdict
        self.df                 =   None
        self.read_meta          =   read_meta
        self.read_map_num       =   int(read_map_num)
        self.read_extra_cols    =   read_extra_cols
        if not self.from_string:
            if filepath is not None:
                self.filepath   =   filepath
                self.filename   =   Path(self.filepath).name
        else:
            self.filepath       =   None
            self.filename       =   Path(filename)
        self.compressed         =   compressed
        self.stat_only          =   stat_only
        self.labels_only        =   labels_only
        self.file_version       =   '0.0'
        self.label_channels     =   0
        self.motifs             =   []
        self.colors             =   []
        self.header             =   ''
        self.default_meta       =   OrderedDict()
        self.extra_comments     =   [
                                        '# Values corresponding to intervals \
                                        (StdDev, HapDelta) refer to the interval \
                                        between current site and next site'
                                    ]
        self.extra_cols         =   []
        self.columns            =   self.default_columns + self.extra_cols
        self.read_om(string=string)
        if self.df is not None:
            if self.mapnum is None:
                if self.mapnum != len(self.df):
                    self.mapnum = len(self.df)
                    warning_msg =   f'Actual map number {len(self.df)}' +\
                                    f'does not match header statement {self.mapnum}.'
                logger.warning(warning_msg)
                self.mapnum = len(self.df)
        else:
            if self.mapnum is None:
                if self.mapnum != len(self.lengths) and not self.stat_only:
                    warning_msg =   f'Actual map number {len(self.lengths)}' +\
                                    f'does not match header statement {self.mapnum}.'
                    logger.warning(warning_msg)
                    self.mapnum =   len(self.lengths)
        logger.info(f'{self.mapnum} maps read.')
        logger.info(f'Complete CMAP file reading.')

    @log_func_run
    def extract_meta(self, string=None):
        """Extract metadata from # comments in a CMAP file"""
        if self.from_string:
            if self.compressed:
                string          =   gzip.decompress(string).decode('utf-8')
            f                   =   io.StringIO(string)
        else:
            if not self.compressed:
                f               =   open(self.filepath, 'r')
            else:
                gz              =   gzip.open(self.filepath, 'r')
                f               =   io.BufferedReader(gz)
        file_ver_line_pattern       =   re.compile(r'[Ff]ile [Vv]ersion:[\t\s]+(\d+\.\d+)')
        label_channels_line_pattern =   re.compile(r'[Ll]abel [Cc]hannels:[\t\s]+(\d+)')
        map_num_line_pattern        =   re.compile(r'[Nn]umber of [Cc]onsensus [Mm]aps:[\t\s]+(\d+)')
        while True:
            line                    =   f.readline()
            if self.compressed and isinstance(line, bytes):
                line = line.decode('utf-8')
            if line[0] != '#':
                logger.debug(f'{line[:10]}...')
                logger.debug('Finished reading metadata.')
                f.close()
                break
            else:
                if self.file_version == '0.0':
                    match_file_ver              =   re.search(file_ver_line_pattern, line[1:].strip())
                if match_file_ver is not None:
                    self.file_version           =   match_file_ver[1]

                    if self.file_version == '0.1':
                        self.default_columns    =   self.default_columns[:9]
                if self.label_channels == 0:
                    match_channels              =   re.search(label_channels_line_pattern, line[1:].strip())
                    if match_channels is not None:
                        self.label_channels     =   int(match_channels[1])
                match_recogn = None
                if len(self.motifs) < self.label_channels:
                    match_recogn                =   re.search('[Nn]ickase [Rr]ecognition [Ss]ite \d'
                                                    + ':[\t\s]+(\S+)\;?(\S+)?', line)
                    if match_recogn is not None:
                        if match_recogn[1]  is not None:
                            self.motifs.append(match_recogn[1])
                        if match_recogn[2]  is not None:
                            self.colors.append(match_recogn[2])
                match_map_num                   =   re.search(map_num_line_pattern, line)
                if match_map_num is not None:
                    self.mapnum                 =   int(match_map_num[1])
                if match_file_ver is None and match_channels is None and \
                    match_recogn is None and match_map_num is None:
                    self.extra_comments.append(line.strip())
                if line[:2] == '#h':
                    self.columns                =   list(map(lambda x: x.strip(),\
                                                    line[3:].split('\t')))
                    if len(self.columns) >= len(self.default_columns):
                        logger.warning('This CMAP file contains non-standard columns.')
                        if not self.read_extra_cols:
                            logger.warning('Dropping info from extra columns.')
                if line[:2] == '#f':
                    self.col_dtype_desc         =   line[2:].split('\t')
        if self.label_channels == 0:
            self.label_channels = 1
        self.default_meta['file_version']       =   '# CMAP File Version:\t' + self.file_version
        logger.info(f'CMAP v{self.file_version}')
        if self.file_version not in self.supported_file_versions:
            logger.warning( f'CMAP gile version {self.file_version} not supported.' +\
                            f'Supported versions are {self.supported_file_versions}')
        self.default_meta['label_channels']     =   str(self.label_channels)
        self.default_meta['mapnum']             =   '# Number of Consensus Maps:\t' + str(self.mapnum)
        for i in range(1,self.label_channels+1):
            if len(self.motifs) < i:
                self.motifs.append('unknown')
            if len(self.colors) < i:
                self.colors.append('unknown')
            self.default_meta['recogn_motif']   =   '# Nickase Recognition Site ' + \
                                                    str(i) + ':\t' + \
                                                    str(self.motifs[i-1]) + ';' + \
                                                    str(self.colors[i-1])
        self.default_meta['h']                  =   '#h\t' + '\t'.join(self.columns)
        self.default_meta['f']                  =   '#f\t' + '\t'.join(self.col_dtype_desc)
        self.header                             =   '\n'.join(self.extra_comments +
                                                    [self.default_meta[x] \
                                                    for x in self.default_meta])
        return self

    @log_func_run
    def read_om(self, string=None):
        """Read CMAP file into CMAP object"""
        if self.read_meta:
            self.extract_meta(string=string)
        if Path(self.filename).suffixes[-1].lower() != '.cmap' and Path(self.filename).suffixes[-2].lower() != '.cmap':
            logger.error(f'{Path(self.filename).suffix.lower()} format is not supported.')
        if self.from_string or string is not None:
            if self.compressed:
                f = io.StringIO(gzip.decompress(string).decode('utf-8'))
        else:
            f = self.filepath

        if self.stat_only:
            df =    pd.read_csv(    f, sep='\t', comment='#',
                                    header=None, index_col=None, usecols=[1,2],
                    ).drop_duplicates()
            df.columns =   ['ContigLength', 'NumSites']
            self.lengths = df['ContigLength'].to_numpy()
            self.mapnum = len(self.lengths)
            self.label_density = float(sum(df['NumSites']))/float(sum(self.lengths)) *100000
            del df
        else:
            self.df     =   pd.read_csv(f, sep='\t', comment='#',
                                header=None, index_col=None, usecols=[0,1,5])\
                            .groupby([0,1]).agg(list).reset_index()
            self.df.columns = ['omid', 'length', 'labels']
            self.df.set_index('omid', inplace=True)
            self.mapnum = len(self.df)
            if not self.labels_only:
                self.omdict     =   {}
                for row in self.df.iterrows():
                    omid = int(row[1]['omid'])
                    self.omdict[omid] = Opmap(  omid=omid,
                                                labels=row[1]['labels'])
            logger.debug(self.df.head())
            self.lengths = self.get_lengths()
        return self

    @log_func_run
    def write_om(self, omdict, filepath=sys.stdout, file_version='0.2'):
        """Write CMAP file output"""
        if filepath != sys.stdout:
            assert Path(filepath).suffix.lower in ['.bnx', '.bnx.gz']
        elif Path(filepath).is_file():
            raise FileExistsError
        # to edit above code that modifies default_columns
        # followed document 'Genome map information black specification' section
        if file_version == '0.2':
            if 'ChimQuality'    not in  self.df.columns:
                self.df['ChimQuality']  =   -1.0
            if 'SegDupL'        not in  self.df.columns:
                self.df['SegDupL']      =   -1.0
            if 'SegDupR'        not in  self.df.columns:
                self.df['SegDupR']      =   -1.0
            if 'FragileL'       not in  self.df.columns:
                self.df['FragileL']     =   0.0
            if 'FragileR'       not in  self.df.columns:
                self.df['FragileR']     =   0.0
            if 'OutlierFrac'    not in  self.df.columns:
                self.df['OutlierFrac']  =   0.0
            if 'ChimNorm'       not in  self.df.columns:
                self.df['ChimNorm']     =   -1.0
            if 'Mask'           not in  self.df.columns:
                self.df['Mask']         =   0
        self.df.columns                 =   self.default_columns
        self.df.to_csv(filepath, comment=self.header, sep='\t')


class RMAP:
    """
    Handle RMAP IO (read file only so far)
    """
    def __init__(self, rmap_filepath, bnx_filepath=None, read_map_num=-1, read_extra_cols=False):
        self.rmap_filepath      =   rmap_filepath
        self.bnx_filepath       =   bnx_filepath
        self.read_map_num       =   read_map_num
        self.read_extra_cols    =   read_extra_cols
        self.default_columns    =   [
                                        'Index', 'MapID', 'Length', 'Color',
                                        'RepeatStart', 'RepeatEnd', 'AvgRepeatLength',
                                        'NRepeatUnits', 'StandardDeviation', 'Confidence',
                                        'RepeatString', 'FPstring', 'FNstring']
        self.col_dtypes         =   [
                                        int, int, float, int, float, float, float, 
                                        int, float, float, str, str, str
                                    ]
        self.read_rmap()

    @log_func_run
    def read_rmap(self):
        """Read RMAP as dataframe"""
        self.df = pd.read_csv(self.rmap_filepath, sep='\t', comment='#', header=None)
        self.df.columns = self.default_columns
        return self

    @log_func_run
    def to_opmaps(self):
        """Extract repeats from RMAP and corresponding BNX as new Opmaps dataset"""
        if self.bnx_filepath is None:
            raise ValueError('Please input BNX file path to extract repeats from maps.')
        bnx             =   BNX(self.bnx_filepath)
        rep_opmaps      =   Opmaps()
        for row in self.df.iterrows():
            rep_om_name =   '_'.join([row['MapID'], row['RepeatStart'], row['RepeatEnd']])
            rep_om_segs =   bnx[row['MapID']][row['RepeatStart']:row['RepeatEnd']]
            rep_opmaps[rep_om_name] = rep_om_segs
        return rep_opmaps


"""
Sequence / in silico digestion
"""

class Seq:
    """
    Sequence class for OM analytic operation
    """
    def __init__(self, seq=None, motif=''):
        """Initialize Seq sequence object"""
        self.seq        =   seq
        self.motif      =   motif

    def rev_comp(self):
        """Get reverse complement of sequence"""
        comp            =   {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        rc              =   ''.join([comp[x] for x in self.seq[::-1]])
        return rc

    def to_opmap(self, fwd=True, rev=True):
        """In silico digestion to create Opmap from sequence; placeholder"""
        motifs          =   {
                                'bsssi': 'CACGAG',
                                'bspqi': 'GCTCTTC',
                                'dle1': 'CTTAAG'
                            }
        if self.motif in motifs:
            self.motif  =   motifs[self.motif]
        f_labels       =   []
        r_labels       =   []
        if fwd:
            f_labels   =   []
        if rev:
            r_labels   =   []
        labels         = sorted(f_labels + r_labels)
        return Opmap(labels=labels)

    def frag_sites(self, thres=200, mode=0):
        """returns 0) both 1) left 2) right 3) mid point of nicking fragile sites; placeholder"""
        sites           =   []
        return sites


class FASTA(object):
    """FASTA sequence dataset object; placeholder"""
    supported_file_formats = ['.fasta', '.fa', 'fas', 'fna']
    def __init__(self, filepath=None):
        self.filepath   =   Path(filepath)
        self.seqs       =   []
        if Path(self.filepath).suffix.lower() not in self.supported_file_formats:
            raise ValueError('Please input a FASTA file.')

    def to_opmaps(self, enzyme):
        """In silico digestion to create Opmaps dataset from FASTA; placeholder"""
        for seq in self.seqs:
            yield Seq(seq).to_opmap()


def calc_n_metrics(num_list, n=50):
    """Calculate N-metrics, by default N50 for input list of numbers."""
    if n>0 and n<100:
        tmp_sum = 0
        total = sum(num_list)
        for l, x in enumerate(sorted(num_list, reverse=True)):
            tmp_sum += x
            if tmp_sum >= total * n/100.0:
                logger.info('N{n}: {result}'.format(n=n, result=x))
                logger.info('L{n}: {result}'.format(n=n, result=l))
                return (x, l+1)
    elif n==100:
        return (min(num_list), n)
    else:
        raise ValueError('N should be a positive number less than 100, generally multiples of 10.')


def labels_to_segments(labels, start=False, use_int=True, omtools=False):
    """Get Opmap segment lengths from label positions"""
    pointer     =   labels[0]
    if start:
        segments = [pointer]
    else:
        segments    =   []
    if omtools:
        for label in labels[1:]:
            segment             =   int(label - 1 - pointer)
            segments.append(segment)
            pointer             +=  segment + 1
    elif not use_int:
        for label in labels[1:]:
            segment             =   round(label - pointer, 1)
            segments.append(segment)
            pointer             +=  segment
    else:
    # to speed up repeat detection
    # for OMTools compatibility, use int without round and -1
    # i.e.,int(label - 1 - pointer)
        for label in labels[1:]:
            segment             =   int(round(label - pointer))
            segments.append(segment)
            pointer             +=  segment
        if segments[-1] ==  19:
            segments[-1] =  20
    return segments


def segments_to_labels(segments, use_int=False):
    """Convert segment length into label positions showing cumulative lengths"""
    labels              =   []
    #labels.append(19)
    pointer             =   0.0
    for segment in segments:
        pointer     += round(segment, 1)
        if  use_int:
            pointer = int(pointer)
        labels.append(pointer)
        
    #labels.append(20)
    return labels


def auto_garbage_collect(pct=80.0):
    """
    auto_garbage_collection - Call the garbage collection if memory used is greater than 80% of total available memory.
                              This is called to deal with an issue in Ray not freeing up used memory.

        pct - Default value of 80%.  Amount of memory in use that triggers the garbage collection call.
        by Michael Wade
        from: https://stackoverflow.com/questions/55749394/how-to-fix-the-constantly-growing-memory-usage-of-ray
    """
    if psutil.virtual_memory().percent >= pct:
        gc.collect()
    return


def filter_bnx(
        job_id=None, lines_to_read=None,
        filepath=None, output_filepath=None, string=None, threads=set_threads(),
        min_length=-1, max_length=-1,
        min_label=-1, max_label=-1,
        min_label_density=-1, max_label_density=-1,
        min_label_snr=-1, max_label_snr=-1,
        min_label_intensity=-1, max_label_intensity=-1,
        min_mol_snr=-1, max_mol_snr=-1,
        min_mol_intensity=-1, max_mol_intensity=-1,
        scan=-1,
    ):
    """On the fly filtering of the input BNX file based on input thresholds"""
    if filepath is not None:
        file_format, compressed =   check_file_extension(filepath)
    if file_format != 'bnx':
        raise ValueError('Input file must by a BNX file')
    bnx = BNX(filepath, compressed=compressed, meta_only=True)
    row0_types = bnx.row0_types
    row0_columns = bnx.row0_columns
    if  compressed:
        gz  =   gzip.open(filepath, 'rb')
        f   =   io.BufferedReader(gz)
    else:
        f = open(filepath)
    output = open(output_filepath, 'w')
    lines = f.readlines()
    mol_in= 0
    logger.info('Start molecule filtering')
    for i, line in enumerate(lines):
        line = line.decode('utf-8')
        if line[0]=='#':
            if 'Number of Molecules' not in line:
                output.write(line.strip()+'\n')
        if line[0] == '0':
            if not compressed:
                qx11 = lines[i+2].strip()
                qx12 = lines[i+3].strip()
            else:
                qx11 = lines[i+2].decode('utf-8').strip()
                qx12 = lines[i+3].decode('utf-8').strip()
            row0_values =   line.split()[2:]
            row0_values =   [assign_dtype(x, (row0_types[2:])[i]) \
                            for i, x in enumerate(row0_values)]
            mol_header  =   dict(zip(row0_columns, row0_values))
            mol_header['LabelDensity'] = mol_header['NumberofLabels']/mol_header['Length']*1e5
            mol_header['LabelSNR'] = np.mean(list(map(float, qx11.split()[1:])))
            mol_header['LabelAvgIntensity'] = np.mean(list(map(float, qx12.split()[1:])))
            if not ((min_length <= mol_header['Length']) and \
                (max_length == -1 or max_length >= mol_header['Length'])):
                continue
            if not ((min_label <= mol_header['NumberofLabels']) and \
                (max_label == -1 or max_label >= mol_header['NumberofLabels'])):
                continue
            if not ((min_label_density <= mol_header['LabelDensity']) and \
                (max_label_density == -1 or max_label_density >= mol_header['LabelDensity'])):
                continue
            if not ((min_mol_snr <= mol_header['SNR']) and \
                (max_mol_snr == -1 or max_mol_snr >= mol_header['SNR'])):
                continue
            if not ((min_mol_intensity <= mol_header['AvgIntensity']) and \
                (max_mol_intensity == -1 or max_mol_intensity >= mol_header['AvgIntensity'])):
                continue
            if not ((min_label_snr <= mol_header['LabelSNR']) and \
                (max_label_snr == -1 or max_label_snr >= mol_header['LabelSNR'])):
                continue
            if not ((min_label_intensity <= mol_header['LabelAvgIntensity']) and \
                (max_label_intensity == -1 or max_label_intensity >= mol_header['LabelAvgIntensity'])):
                continue
            if not (scan == -1 or scan == mol_header['ScanNumber']):
                continue
            output.write('\n'.join([line.strip(), lines[i+1].decode('utf-8').strip(), qx11, qx12])+'\n')
            mol_in += 1
        if i%10000 == 0:
            auto_garbage_collect()
        if i%1e5 == 0:
            logger.info(f'{i:.0f} file lines read')
            logger.info(f'{mol_in} molecules filtered in')
    f.close()
    output.close()
    logger.info(f'Finished BNX filtering. {mol_in} molecules remain.')
    return 


def filter_bnx_chunk_wrapper(
        job_id, lines_to_read, row0_types, row0_columns,
        filepath, output_filepath, string, threads,
        min_length, max_length,
        min_label_density, max_label_density,
        min_label_snr, max_label_snr,
        min_label_intensity, max_label_intensity,
        min_mol_snr, max_mol_snr,
        min_mol_intensity, max_mol_intensity, q,
):
    auto_garbage_collect()
    return q.put(filter_bnx_chunk(
                    job_id, lines_to_read, row0_types, row0_columns,
                    filepath, output_filepath, string, threads,
                    min_length, max_length,
                    min_label_density, max_label_density,
                    min_label_snr, max_label_snr,
                    min_label_intensity, max_label_intensity,
                    min_mol_snr, max_mol_snr,
                    min_mol_intensity, max_mol_intensity,
    ))


def filter_bnx_chunk(
        job_id, lines_to_read, row0_types, row0_columns,
        filepath, output_filepath, string, threads,
        min_length, max_length,
        min_label_density, max_label_density,
        min_label_snr, max_label_snr,
        min_label_intensity, max_label_intensity,
        min_mol_snr, max_mol_snr,
        min_mol_intensity, max_mol_intensity,
):
    """Filter BNX file"""
    tmp_bnxstat_list    =   []
    for i, line in enumerate(lines_to_read):
        if i!=0 and not i%10000:
            logger.debug(f'{line[:10]}...')
            logger.info(f'{i} molecules read in job {job_id}')
        split_line = line.split()
        omid = int(split_line[0])
        if split_line[1] == '0':
            tmp_bnxstat_list.append([])
            row0_values =   ['1'] + split_line[2:]
            row0_values =   [assign_dtype(x, (row0_types)[i]) \
                            for i, x in enumerate(row0_values)]
            mol_header  =   dict(zip(row0_columns, row0_values))
            mol_header['MoleculeId'] = omid
            chip = mol_header['ChipId'].strip('chips,')
            tmp_bnxstat_list[-1] += [
                mol_header['MoleculeId'],
                mol_header['Length'],
                mol_header['NumberofLabels'],
                mol_header['SNR'],
                mol_header['AvgIntensity'],
                mol_header['Flowcell'],
                mol_header['ScanNumber'],
                -1,
                -1
            ]
        elif split_line[1] == '1':
            try:
                labels = list(map(float, split_line[2:]))
                tmp_bnxstat_list[-1][9] = labels
            except IndexError:
                tmp_bnxstat_list.append([omid]+[-1]*8)
                tmp_bnxstat_list[-1] += labels
        elif split_line[1][0] == 'Q':
            arr = np.array(list(map(float, split_line[2:])))
            if split_line[1] == 'QX11':
                if len(split_line) == 2:
                    qx11_mol_mean = 0
                else:
                    qx11_mol_mean = numba_mean(arr)
                try:
                    tmp_bnxstat_list[-1][7] = qx11_mol_mean
                except IndexError:
                    tmp_bnxstat_list.append([omid]+[-1]*6)
                    tmp_bnxstat_list[-1] += [qx11_mol_mean, -1]
            if split_line[1] == 'QX12':
                if len(split_line) == 2:
                    qx12_mol_mean = 0
                else:
                    qx12_mol_mean = numba_mean(arr)
                try:
                    tmp_bnxstat_list[-1][8] = qx12_mol_mean
                except IndexError:
                    tmp_bnxstat_list.append([omid]+[-1]*7)
                    tmp_bnxstat_list[-1] += [qx11_mol_mean]
        elif split_line[1][0] not in ['#', '0', '1', 'Q']:
                logger.warning(f'Invalid line detected: {line[:20]}...')
                logger.warning(f'A line in BNX file should start with #, 0, 1 or QX, but not {split_line[1]}')
    auto_garbage_collect()
    return tmp_bnxstat_list

@jit(cache=True, nopython=True, parallel=False)
def numba_mean(a):
    return np.mean(a)


@jit(cache=True, nopython=True, parallel=False)
def numba_std(a):
    return np.std(a)