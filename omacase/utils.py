#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220228


"""
Utility functions for OMCASE
By Claire Chung
"""


import logging
from datetime import datetime
import inspect
from pathlib import Path

from deepmerge import Merger
import jinja2

from omacase.param_parser import ALIGN_SUPPORTED_FILE_FORMATS, ANNO_SUPPORTED_FILE_FORMATS
from omacase.param_parser import CONFIG_SUPPORTED_FILE_FORMATS, OPMAPS_SUPPORTED_FILE_FORMATS
from omacase.param_parser import DEFAULT_COLORSCHEME

logger = logging.getLogger(__name__)

LOG_FORMAT = "%(asctime)s [%(threadName)-12.12s] [%(module)-10.10s] [%(levelname)-5.5s] %(message)s"

def banner_art():
    """Prints banner art"""
    print(" ,^\_,^\ ,--.   ,--.  ,^-^.   ,-----.   ,^-^.    ,---.   ,-----.")
    print("'  • •  '|   `.'   | /  ▽  \ '- .--.'  /  ▽  \  '-.--.' |- .---'")
    print("|-  ▽  -||- |'.'| -||- .-. -||- |     |- .-. -| \-\--., |- `--, ,-,")
    print("'- '-' -'|- |   | -||- | | -|'- '__.-,|- | | -|'-'.__' '|- `--.' ,'")
    print(" `-----' `--'   `--'`--' `--' `-----' `--' `--' `-----' `------'")
    return


def check_file_extension(filepath, datatype='opmaps'):
    """Check file extension."""
    path_suffixes = Path(str(filepath.lower())).suffixes
    logger.debug(f'{path_suffixes}')
    if path_suffixes[-1] == '.gz':
        file_format, compressed =   path_suffixes[-2].strip('\.'), True
        logger.info('Input file is GZ compressed')
    else:
        file_format, compressed =   path_suffixes[-1].strip('\.'), False
    if datatype == 'opmaps':
        supported = OPMAPS_SUPPORTED_FILE_FORMATS
    elif datatype == 'annotation':
        supported = ANNO_SUPPORTED_FILE_FORMATS
    elif datatype == 'alignment':
        supported = ALIGN_SUPPORTED_FILE_FORMATS
    if file_format not in supported:
        err =   f'Input file: {filepath}\n' + \
                f'Invalid file format: {file_format}\n\
                Supported file formats: {supported}'
        logger.error(err)
        raise   ValueError(err)
    return  file_format, compressed


def get_file_segment_start_byte(file_bytes, bytes_per_chunk):
    return [b - 1 for b in range(0, file_bytes, bytes_per_chunk)]

"""
def read_file_segment(f, start_byte, bytes_per_chunk):
    f.seek(start_byte)
    bytes_read = f.read(bytes_per_chunk - start_byte + 1)
    logger.debug('to pickle')
    return pickle.dumps(bytes_read)
"""

def log_func_run(func):
    def wrapper(*args, **kwargs):
        module = Path(inspect.stack()[1].filename).stem
        logger.info("Start running {}.{}".format(module, func.__name__))
        func(*args, **kwargs)
        logger.info("Complete running {}.{}".format(module, func.__name__))
    return wrapper


def get_stack(func):
    def wrapper(*args, **kwargs):
        stack = inspect.stack()
        logger.debug("Stack of {}: {}".format(func.__name__, stack))
        func(*args, **kwargs)
    return wrapper


def assign_dtype(var, dtype):
    var = dtype(var)
    return var


def custom_round(x, base=200):
    return base * round(x/base)


def get_layout(layout_dict=None, **kwargs):
    if 'colors' in kwargs:
        color_scheme = get_color_scheme(kwargs['colors'])
    else:
        color_scheme = get_color_scheme(DEFAULT_COLORSCHEME)
    layout_default = {
        'yaxis':            {'visible':         False,
                             'showticklabels':  False,
                             'range':           [0,6]
                            },
        'hoverlabel':       {'font_size': 16},
        'hoverdistance':    100,
        'margin':           {'l': 20, 'r': 20, 't':20, 'b':20},
        'paper_bgcolor':    color_scheme['paper_bg']
    }
        
    merger  =   Merger([
                    (dict, ["merge"])
                    ],
                    ["override"], ["override"]
                )
    layout = merger.merge(layout_default, kwargs)
    if layout_dict:
        layout = merger.merge(layout, layout_dict)
    return layout

def get_color_scheme(colors):
    color_schemes = {
        # dark background, light foreground
        'yoyo1gfp':     {   'fillcolor':        'Navy',
                            'line_color':       'MediumSpringGreen',
                            'highlight_color':  'Yellow'},
        'yoyo1rfp':     {   'fillcolor':        'Navy',
                            'line_color':       'Salmon',
                            'highlight_color':  'Yellow'},
        'subdue':       {   'fillcolor':        'DarkSlateGray',
                            'line_color':       'Gainsboro',
                            'paper_bg':         'LightSteelBlue',
                            'highlight_color':  'Gold'},
        'test0':         {   'fillcolor':       '#1A2C56', # dark blue
                            'line_color':       '#00CC96', # Carribean Green
                            'paper_bg':         'LightSteelBlue',
                            'highlight_color':  'Gold'},
        'test1':         {   'fillcolor':       '#1A2C56',
                            'line_color':       'Peachpuff',
                            'paper_bg':         'LightSteelBlue',
                            'highlight_color':  'Powderblue'},
        'grape':       {   'fillcolor':        'DarkSlateBlue',
                            'line_color':       'MistyRose',
                            'paper_bg':         'LightSteelBlue',
                            'highlight_color':  'Plum'},
        'mint':         {   'fillcolor':        'DarkSlateGray',
                            'line_color':       'LightCyan',
                            'paper_bg':         'LightSteelBlue',
                            'highlight_color':  'MediumAquamarine'},
        'buttercup':    {   'fillcolor':        'DarkSlateGray',
                            'line_color':       'HoneyDew',
                            'paper_bg':         'LightSteelBlue',
                            'highlight_color':  'Gold'},
        'purple':       {   'fillcolor':        'Purple',
                            'line_color':       'Orchid',
                            'highlight_color':  'Yellow'},
        'starrynight':  {   'fillcolor':        'MidnightBlue',
                            'line_color':       'Gold',
                            'highlight_color':  'SkyBlue'},

        # light background, dark foreground
        'omtools':      {   'fillcolor':    'Yellow',
                            'line_color':   'Black'}
    }
    return color_schemes[colors]


def render_template_without_request(template_name, **template_vars):
    """
    Usage is the same as flask.render_template:

    render_template_without_request('my_template.html', var1='foo', var2='bar')
    https://stackoverflow.com/questions/17206728/attributeerror-nonetype-object-has-no-attribute-app
    Credits to @chucksmash
    """
    env = jinja2.Environment(
        loader=jinja2.PackageLoader(package_name='omacase', package_path='templates')
    )
    template = env.get_template(template_name)
    return template.render(**template_vars)


def timestamp(format="%Y%m%d-%H%M%S"):
    return f'{datetime.now().strftime(format)}'

def to_list(*args):
    return([arg for arg in args])

def easter_cat():
    """Returns easter egg cat, meow~"""
    cat = "ฅ(ﾐ・ﻌ・ﾐ)ฅ"
    return cat
