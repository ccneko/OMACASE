#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220219


"""
OM data web viewer
By Claire Chung
"""


import logging

from deepmerge import Merger
from  plotly.subplots import make_subplots
from dash import dcc, html, Output, Input, State#, Download
import dash_bootstrap_components as dbc

from omacase.opmap import Opmap, Opmaps, set_threads
from omacase.alignment import OMAlignments
from omacase.annotations import get_anno, Annotations, GTF
from omacase.dataviz import *
from omacase.param_parser import MAPS_PER_PAGE, MAX_X, DEFAULT_COLORSCHEME


logger = logging.getLogger(__name__)


class DisplayAnnotations:
    """Display annotations"""
    # intervaltree of intervals & annotations (BED/GTF/GFF)
    def __init__(   self,
                    args=None,
                    r_opmaps=None,
                    q_opmaps=None,
                    filepath=None,
                    annotations: Annotations=None,
                    alignments=None,
                    app=None
        ):
        filepath = '/home/claire/R/x86_64-pc-linux-gnu-library/4.0/GenomicFeatures/extdata/GTF_files/Aedes_aegypti.partial.gtf'
        if annotations is not None:
            self.annotations    =   annotations
        elif annotations is None and filepath is not None:
            self.annotations = get_anno(filepath=filepath)
        elif not isinstance(annotations, Annotations):
            raise ValueError(
                f'Input (type {type(annotations)})' + \
                'is not Annotations dataset object.')
        else:
            logger.warning('No annotation was imported. Please import a dataset to proceed.')
        self.r_opmaps       =   r_opmaps
        self.q_opmaps       =   q_opmaps
        self.alignments     =   alignments
        self.layout         =   self.display_annotations()

    def display_annotations(self):
        return dcc.Graph(
            figure  =   draw_anno(
                            annotations=self.annotations,
                            alignments=self.alignments,
                            r_opmaps = self.r_opmaps,
                            q_opmaps = self.q_opmaps,
                        )) # GTF only now


class DisplayOMAlignments:
    """Display OM Alignments"""
    def __init__(self, r_opmaps: Opmaps, q_opmaps: Opmaps, oma: OMAlignments):
        raise NotImplemented
        self.r_opmaps       =   r_opmaps
        self.q_opmaps       =   q_opmaps
        self.oma            =   oma
        self._app        =   dash.Dash(__name__)
        self._app_layout()
        self.display_alignments()

    def display_alignments( self, ref_id: int = 1, 
                            start_coord: int = 0, end_coord: int = -1,
                            start_label: int = None, end_label: int = None
                        ):
        raise NotImplemented
        draw_om(r_opmaps[ref_id], start_coord=0, end_coord=-1)
        # ref = draw_om(opmap,start_label=, end_label=100)
        return self