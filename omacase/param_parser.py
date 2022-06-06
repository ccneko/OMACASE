#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20210507


"""
Parse options from command-line, JSON or config.ini
"""


__author__ = 'Claire Chung'


import argparse
import json
from pathlib import Path
import platform
import multiprocessing as mp

import toml

OPMAPS_SUPPORTED_FILE_FORMATS   =   ['bnx', 'cmap']
ALIGN_SUPPORTED_FILE_FORMATS    =   ['xmap', 'oma']
ANNO_SUPPORTED_FILE_FORMATS     =   ['bed', 'gtf']
CONFIG_SUPPORTED_FILE_FORMATS   =   ['json', 'toml']

## Data reader settings
THREADS             =   mp.cpu_count()

## Web app settings
PORT                =   8050
MAX_DATASETS        =   5

## Viewer settings
MAPS_PER_PAGE       =   10
SORT_BY             =   'id'
REPEAT_ONLY         =   False
MAX_X               =   1e6
REVERSE             =   False
DEFAULT_COLORSCHEME =   'buttercup'

## Repeat calling parameters
MAX_DIFF            =   800
MAX_DIFF_RATIO      =   0.1
MAX_AVG_DIFF        =   1000
MAX_AVG_DIFF_RATIO  =   0.08
K                   =   13
MIN_COPY_NUM        =   5
MIN_SEG_LEN         =   500
MIN_CASETTE_LENGTH  =   2000
MAX_MISSED_LABEL    =   1
RESOLUTION          =   250

## Repeat simulator parameters
## check repeat_simulator.py for more description
RANDOM_SEED = 22
SAMPLE_N = 10
REP_COUNT_MAX = 3
REP_COUNT_MIN = 0
REP_COUNT_SKEW = 4
REP_COUNT_LOC = 0
REP_COUNT_SCALE = 1
K_MIN = 3
K_MAX = 13
K_SKEW = 4
K_LOC = 3
K_SCALE = 3
CN_MIN = 5
CN_MAX = 50
CN_SKEW = 4
CN_LOC = 5
CN_SCALE = 3
LEN_ERR_MAX = 0.03
REPEAT_UNIT_LENGTH_MIN = 3000
REPEAT_UNIT_LENGTH_MAX = 500000
SEG_LENGTH_MIN = 500
SEG_LENGTH_MAX = 50000
TOTAL_REPEAT_LENGTH_MAX = 1000000

## Uploader settings
machine_platform = platform.system()
if machine_platform in ['Linux', 'Darwin']:
    UPLOAD_FOLDER_ROOT  = '/tmp/omacase-uploads/'
elif machine_platform in ['Windows']:
    UPLOAD_FOLDER_ROOT = "C:\\tmp\omacase-uploads\\"
## Log settings
VERBOSE             =   1

def config_file_from_args(args, filepath='config.toml'):
    filepath = Path(filepath)
    file_extension = filepath.suffix[1:]
    if file_extension not in CONFIG_SUPPORTED_FILE_FORMATS:
        raise ValueError(   f'Invalid file type {file_extension}. Supported formats are \
                            {list(map(lambda x: x.upper, CONFIG_SUPPORTED_FILE_FORMATS))}.')
    if file_extension == 'json':
        with open(filepath, 'w') as f:
            json.dump(args, f, indent=4)
            return json.dumps(args, f, indent=4)
    elif file_extension == 'toml':
        return toml.dump(args, open(filepath, 'w'))


def parse_opts_from_file(filepath=''):
    """Get options from config files (JSON or TOML)"""
    filepath = Path(filepath)
    file_extension = filepath.suffix[1:]
    if file_extension not in CONFIG_SUPPORTED_FILE_FORMATS:
        raise ValueError(   f'Invalid file type {file_extension}. Supported formats are \
                            {list(map(lambda x: x.upper, CONFIG_SUPPORTED_FILE_FORMATS))}.')
    if file_extension == 'json':
        with open(filepath, 'rb') as f:
            args = json.load(f)
    elif file_extension == 'toml':
        args = toml.load(filepath)
    return args


#globals().update(config_file_from_args('config/default.toml'))


def parse_opts_from_cmd(args=[]):
    """Get options from command line"""
    parser = argparse.ArgumentParser(   prog='OMACASE',
                                        description='Optical mapping toolkit')
    parser.add_argument('-c', '--create_config', default='config.toml', nargs='?',
                        help='Dry run to create configuration file from arguments')
    either          =   parser.add_mutually_exclusive_group()
    either.add_argument('-d', '--input_directory', metavar='INPUT_DIR_PATH',
                        type=str, default=None,
                        help='Input folder path')
    either.add_argument('-i', '--input_file', metavar='INPUT_FILE_PATH',
                        type=str, default=None, nargs='?',
                        help='Input file path')
    parser.add_argument('mode', type=str, default='web',
                        choices = [
                                    'web', 'w',
                                    'view', 'v',
                                    'report', 'r',
                                    'text-report', 't',
                                    'repeat-bed', 'b',
                                    'repeat-sim', 's',
                                    '',
                        ],
                        metavar='MODE',
                        help='Select OMACASE mode: web, view, report, ' + \
                            'text-report, repeat-bed, repeat-sim')
    parser.add_argument('-e', '--excluded_regions_bed', default=None, nargs='?',
                        help='regions to exclude from simulation specified in BED file')
    parser.add_argument('-f', default=None, nargs='?', metavar='IGNORE_THIS',
                        type=str, help='Placeholder for running on Jupyter')
    parser.add_argument('-F', '--report_output_format', default=None, nargs='*',
                        choices = ['pdf', 'tsv'],
                        help='Output formats available: pdf (provides plots), tsv (parsable)')
    parser.add_argument('-n', '--read_map_num', type=int, default=-1,
                        metavar='N', help='First N maps to read')
    parser.add_argument('-N', '--maps_per_page', default=MAPS_PER_PAGE, nargs='?',
                        type=int, metavar='N', help='Number of maps to show on a page')
    parser.add_argument('-o', '--output_prefix', default=None, nargs='?', type=str,
                        help='Output report prefix')
    parser.add_argument('-p', '--port', default=PORT, nargs='?',
                        type=int, metavar='N', help='Binding port for web socket')
    parser.add_argument('-r', '--repeat_only', const=True, nargs='?', default=False,
                        help='Show repeat-containing maps only.')
    parser.add_argument('-R', '--rmap', default=None, nargs='?',
                        help='Input Bionano RMAP format of repeat data to include \
                        visualizations in report.')
    parser.add_argument('-s', '--seed', default=RANDOM_SEED, nargs='?',
                        type=int, help='Random seed for simulation')
    parser.add_argument('-t', '--threads', default=THREADS, nargs='?',
                        type=int, metavar='N', help='Max number of threads to run')
    parser.add_argument('-T', '--timestamp', const=True, nargs='?', default=False,
                        help='Add timestamp to report output file name prefix')
    parser.add_argument('-v', '--verbose', const=True, nargs='?', default=False,
                        help='Switch to debug log level')
    parser.add_argument('--read_meta', const=True, nargs='?')
    parser.add_argument('--read_extra_cols', const=False, nargs='?')
    args = parser.parse_args()
    return args
