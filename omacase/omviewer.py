#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220419


"""
OM data web viewer
By Claire Chung
"""


import logging
import itertools
import math

import plotly.express as px
import plotly.graph_objects as go
from  plotly.subplots import make_subplots
from dash import dcc, html, callback_context, Input, Output, State
import dash_bootstrap_components as dbc

from omacase.app import app
from omacase.annotations import Annotations
from omacase.layouts import QuantitySetter, dummy_config_html
from omacase.opmap import Opmaps, get_opmaps
from omacase.param_parser import *
from omacase.repeat_caller import find_tandem_repeats_from_opmaps_ivtree, get_repeat_score
from omacase.utils import get_color_scheme, get_layout

logger = logging.getLogger(__name__)


class DisplayMaps:
    """Display maps in an OM dataset"""
    def __init__(   self,
                    args                =   None,
                    opmaps: Opmaps      =   None,
                    app                 =   None,
                    threads             =   THREADS,
                    maps_per_page       =   MAPS_PER_PAGE,
                    repeat_only         =   False,
                    sort_by             =   SORT_BY,
                    max_x               =   MAX_X,
                    reverse             =   REVERSE
                ):
        logger.debug(f'outside: {opmaps}')
        self.opmaps = None
        if app is not None:
            self._app = app
        if opmaps is not None:
            self.opmaps         =   opmaps
            logger.debug(f'opmaps is not None: {self.opmaps}')
        if args is not None:
            self.threads        =   args.threads
            THREADS             =   args.threads
            self.maps_per_page  =   args.maps_per_page
            MAPS_PER_PAGE       =   args.maps_per_page
            self.repeat_only    =   args.repeat_only
            if not self.opmaps and not opmaps:
                self.opmaps     =   get_opmaps(args=args, labels_only=True)
        else:
            self.threads        =   threads
            self.maps_per_page  =   maps_per_page
            self.repeat_only    =   repeat_only
        if self.opmaps:
            if not isinstance(self.opmaps, Opmaps):
                raise ValueError(f'Input (type {type(self.opmaps)}) is not Opmaps OM dataset object.')
            self.omids_to_draw  =   self.opmaps.df.index
            omids_to_draw       =   self.omids_to_draw[:self.maps_per_page]
        else:
            logger.warning('No maps input.')
        if self.maps_per_page < 1:
            logger.warning(f'Invalid maps_per_page input {self.maps_per_page}. Revert to default.')
            self.maps_per_page  =   MAPS_PER_PAGE
        self.sort_by            =   sort_by
        self.max_x              =   max_x
        self.reverse            =   reverse
        self.total_page         =   self.get_total_page()
        logger.debug(f'init total page: {self.total_page}')
        self.curr_page          =   1
        self.layout             =   []
        self.repeat_detected    =   False
        self.view_all           =   True
        if self.repeat_only and not self.repeat_detected:
            self.view_all       =   False
            self.call_repeats()
            if len(self.omids_to_draw) > 0:
                self.set_repeat_view()
            self.repeat_detected = True
        if self.repeat_only and self.repeat_detected:
            self.total_page     =   self.get_total_page()
        elif not self.repeat_only and len(omids_to_draw)>0:
            max_map_len = round(max(self.opmaps.df.loc[omids_to_draw, 'labels']\
                                .apply(lambda x:x[-1])), -4)
            max_len = min(max_map_len, self.max_x)
            logger.debug(f'max_len: {max_len}')
            self.view_region    =   [-1e5, max_len]
        logger.debug(f'DisplayMaps.__init__->self.omids_to_draw: {self.omids_to_draw}')
        self.layout             =   self.get_app_layout()
        if self._app is not None and hasattr(self, "callbacks"):
            self.callbacks(self._app)
    
    def callbacks(self, _app):
        @_app.callback(
                    Output('page', 'value'),                    #1
                    Output('total_page', 'value'),              #2
                    Output('by_id', 'active'),                  #3
                    Output('by_length', 'active'),              #4
                    Output('by_repeat_score', 'active'),        #5
                    Output('reverse', 'active'),                #6
                    Output('repeat_only', 'active'),            #7
                    Output('custom_omids', 'value'),            #8
                    Output('omviewer', 'children'),             #9
                    Output('view_region_min', 'value'),         #10
                    Output('view_region_max', 'value'),         #11
                    Output('view_all', 'active'),               #12

                    Input('custom_omids', 'n_submit'),          #1
                    Input('custom_omids_submit', 'n_clicks'),   #2
                    Input('prev_page', 'n_clicks'),             #3
                    Input('next_page', 'n_clicks'),             #4
                    Input('page', 'n_submit'),                  #5
                    Input('goto_page', 'n_clicks'),             #6
                    Input('maps_per_page', 'value'),            #7
                    Input('view_region_min', 'value'),          #8
                    Input('view_region_max', 'value'),          #9
                    Input('view_region_min', 'n_submit'),       #10
                    Input('view_region_max', 'n_submit'),       #11
                    Input('by_id', 'n_clicks'),                 #12
                    Input('by_length', 'n_clicks'),             #13
                    Input('by_repeat_score', 'n_clicks'),       #14
                    Input('reverse', 'n_clicks'),               #15
                    Input('repeat_only', 'n_clicks'),           #16
                    Input('view_all', 'n_clicks'),              #17
                    
                    State('page', 'value'),                     #1
                    State('by_id', 'active'),                   #2
                    State('by_length', 'active'),               #3
                    State('by_repeat_score', 'active'),         #4
                    State('reverse', 'active'),                 #5
                    State('repeat_only', 'active'),             #6
                    State('custom_omids', 'value'),             #7
                    prevent_initial_call=True,
        )
        def update_omviewer(custom_omids_input_n, custom_omids_submit_n,
                            prev_n, next_n, page_n, goto_n,
                            maps_per_page, view_region_min, view_region_max,
                            view_region_min_n, view_region_max_n,
                            sort_by_id_n, sort_by_len_n, sort_by_repeat_score_n,
                            reverse_n, repeat_only_n, view_all_n,
                            curr_page, sort_by_id, sort_by_len, sort_by_repeat_score, reverse,
                            repeat_only, custom_omids,
            ):
            changed_id = [p['prop_id'] for p in callback_context.triggered][0]
            logger.info(f'Last changed element ID: {changed_id}')
            logger.debug(f'self.omids_to_draw: {self.omids_to_draw}')
            if sort_by_id:
                self.sort_by = 'id'
            elif sort_by_len:
                self.sort_by = 'length'
            elif sort_by_repeat_score:
                self.sort_by = 'repeat_score'
            if 'prev_page' in changed_id or \
                'next_page' in changed_id or \
                'goto_page' in changed_id:
                if 'prev_page' in changed_id:
                    turn        =   -1
                    logger.info('Turn to the previous page')
                    to_page     =   max(1, int(curr_page) + turn)
                    to_page     =   min(to_page, self.total_page)
                elif 'next_page' in changed_id:
                    turn        =   1
                    logger.info('Turn to the next page')
                    to_page     =   max(1, int(curr_page) + turn)
                    to_page     =   min(to_page, self.total_page)
                elif 'page.n_submit' in changed_id or 'goto_page' in changed_id:
                    to_page     =   curr_page
                self.curr_page  =   to_page
                if not self.repeat_only:
                    self.view_region = [0, self.max_x]
                else:
                    self.set_repeat_view()
                logger.info(f'Go to page {self.curr_page}')
            elif 'view_region_min' in changed_id or 'view_region_max' in changed_id:
                self.omids_to_draw  =   custom_omids
                if view_region_min != 'None' and view_region_min is not None:
                    logger.debug(f'view_region_min: {view_region_min}')
                    self.view_region = [int(view_region_min), self.view_region[1]]
                if view_region_max != 'None' and view_region_max is not None:
                    self.view_region = [self.view_region[0], int(view_region_max)]
                    logger.debug(f'view_region_max: {view_region_max}')
            elif 'maps_per_page' in changed_id:
                self.maps_per_page  =   maps_per_page
                self.sort_by        =   'id'
                sort_by_id          =   True
                sort_by_len         =   False
                sort_by_repeat_score =  False
                self.reverse        =   False
                self.curr_page      =   1
            elif 'by_id' in changed_id or \
                'by_length' in changed_id or \
                'by_repeat_score' in changed_id:
                if 'by_id' in changed_id:
                    self.sort_by    =   'id'
                    sort_by_id      =   True
                    sort_by_len     =   False
                    sort_by_repeat_score = False
                    self.reverse    =   False
                elif 'by_length' in changed_id:
                    self.sort_by    =   'length'
                    sort_by_id      =   False
                    sort_by_len     =   True
                    sort_by_repeat_score = False
                    self.reverse    =   False
                elif 'by_repeat_score' in changed_id:
                    self.sort_by    =   'repeat_score'
                    sort_by_id      =   False
                    sort_by_len     =   False
                    sort_by_repeat_score = True
                    self.reverse    =   True
                self.curr_page      =   1
                if not self.repeat_only:
                    self.view_region =  [0, self.max_x]
            elif 'reverse' in changed_id:
                self.reverse        =   not self.reverse
                if sort_by_id:
                    self.sort_by    =   'id'
                    sort_by_len     =   False
                elif sort_by_len:
                    self.sort_by    =   'length'
                    sort_by_id      =   False
                self.curr_page      =   1
                if not self.repeat_only:
                    self.view_region =  [0, self.max_x]
            elif 'repeat_only' in changed_id:
                self.repeat_only    =   not self.repeat_only
                if self.repeat_only:
                    if not self.repeat_detected:
                        self.call_repeats()
                    self.set_repeat_view()
                    self.view_all   =   False
                    logger.debug(f'Display repeat-containing maps only: {self.omids_to_draw}')
                    custom_omids    =   ','.join(map(str, self.omids_to_draw))
                    
                else:
                    self.view_all   =   True
                    self.omids_to_draw = self.opmaps.df.index
                    self.view_region =  [0, self.max_x]
                self.sort_by        =   'id'
                sort_by_id          =   True
                sort_by_len         =   False
                self.reverse        =   False
                self.curr_page      =   1
            elif 'custom_omids' in changed_id:
                self.view_all       =   False
                logger.info(f'Custom Opmap IDs submitted: {custom_omids}')
                logger.debug(f'type(custom_id): {type(custom_omids)}')
                self.omids_to_draw  =   custom_omids
                self.curr_page      =   1
                self.view_region    =   [0, self.max_x]
            elif 'view_all' in changed_id:
                self.view_all       =   True
                self.omids_to_draw  =   self.opmaps.df.index
                self.sort_by        =   'id'
                sort_by_id          =   True
                sort_by_len         =   False
                self.reverse        =   False
                self.repeat_only    =   False
                custom_omids        =   'all'
                self.curr_page      =   1
                self.view_region    =   [0, self.max_x]
                self.maps_per_page  =   MAPS_PER_PAGE
            omids_on_page = self.get_target_omids()
            maps_drawn = [self.draw_maps(omids_on_page)]
            self.total_page = self.get_total_page()
            logger.debug(f'self.omids_to_draw: {self.omids_to_draw}')
            logger.debug(f'omids_on_page: {omids_on_page}')
            logger.debug(f'total_page: {self.total_page}')
            logger.debug(f'sort by id: {sort_by_id}')
            logger.debug(f'sort by length: {sort_by_len}')
            return  str(self.curr_page), str(self.total_page), \
                    sort_by_id, sort_by_len, sort_by_repeat_score, self.reverse, self.repeat_only, \
                    custom_omids, maps_drawn, self.view_region[0], self.view_region[1], \
                    self.view_all

    def get_app_layout(self):
        maps_drawn = self.draw_maps(self.get_target_omids())
        layout = html.Div([
            dcc.Location(id='url_viewer', href='/map-viewer/', refresh=False),
            dcc.Loading(
                id          =   'page_loading',
                type        =   'default',
                children    =   [
                    html.Div(
                        id='mapviewer_content', 
                        children=[
                            html.Div(
                                id='filename',
                                className='row',
                                children=html.H3(f'File: {self.opmaps.filename}')
                            ),
                            html.Div(
                                id='viewer_control', 
                                className='row', children=[
                                    html.Div(
                                        className='row',
                                        children=[
                                            html.Div(
                                                'Page control', 
                                                className='col-md-2 col-xs-12 control-labels'),
                                            html.Div(
                                                'Current page', 
                                                className='col-md-4 col-xs-12 control-labels'),
                                            html.Div(
                                                'Maps per page',
                                                className='col-md-2 col-xs-12 control-labels'),
                                            html.Div(
                                                'View region',
                                                className='col-md-4 col-xs-12 control-labels'),
                                        ]
                                    ),
                                    html.Div(className='row', children=[
                                        html.Div(
                                            id = 'pagination',
                                            children = [
                                                html.Div(
                                                    className='col-md-2 col-xs-12',
                                                    children=[
                                                    self.prev_page_button(),
                                                    self.next_page_button(),
                                                    ],
                                                ),
                                                html.Div(
                                                    className='col-md-4 col-xs-12',
                                                    children=[
                                                        self.goto_page_input(),
                                                        html.Div('/', className='li-separator'),
                                                        self.total_page_input(self.total_page),
                                                        self.goto_page_button(),
                                                    ],
                                                ),
                                                html.Div(
                                                    className='col-md-2 col-xs-12',
                                                    children=[
                                                        self.maps_per_page_input(),
                                                    ],
                                                ),
                                                html.Div(
                                                    className='col-md-4 col-xs-12',
                                                    children=[
                                                        self.view_region_min_input(),
                                                        html.Div('-', className='li-separator'),
                                                        self.view_region_max_input(),
                                                    ],
                                                ),
                                            ],
                                        ),
                                    ]),
                                    html.Div(
                                        className='row',
                                        children=[
                                                html.Div(
                                                    'Filter by Map IDs',
                                                    className='col-xl-6 col-md-5 col-xs-12 control-labels',
                                                ),
                                                html.Div(
                                                    'Other filtering',
                                                    className='col-xl-2 col-md-3 col-xs-12 control-labels',
                                                ),
                                                html.Div(
                                                    'Sort by',
                                                    className='col-xl-3 col-md-3 col-xs-12 control-labels',
                                                ),
                                                html.Div(
                                                    'Order',
                                                    className='col-xl-1 col-md-1 col-xs-12 control-labels',
                                                ),
                                        ],
                                    ),
                                    html.Div(
                                        className='row',
                                        children=[
                                            html.Div(
                                                id = 'filter_sort_omids',
                                                children = [
                                                    html.Div(
                                                        className='col-xl-6 col-md-5 col-xs-12',
                                                        children=[
                                                            self.custom_id_input(),
                                                            self.custom_id_submit(),
                                                        ],
                                                    ),
                                                    html.Div(
                                                        className='col-xl-2 col-md-3 col-xs-12',
                                                        children=[
                                                            self.repeat_only_button(),
                                                            self.view_all_button(),
                                                        ],
                                                    ),
                                                    html.Div(
                                                        className='col-xl-3 col-md-3 col-xs-12',
                                                        children=[
                                                            self.sort_by_button_group(),
                                                        ],
                                                    ),
                                                    html.Div(
                                                        className='col-xl-1 col-md-1 col-xs-12',
                                                        children=[
                                                            self.sort_reverse_button(),
                                                        ],
                                                    ),
                                                ]
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                            html.Div(id='omviewer', children = maps_drawn),
                            dummy_config_html,
                        ]
                    ),
                ]
            )
        ])
        return layout

    def call_repeats(self):
        """Call repeats and return results to viewer"""
        self.repeat_dict    =   find_tandem_repeats_from_opmaps_ivtree(
                                    opmaps = self.opmaps,
                                    threads = self.threads,
                                )
        logger.info(f'Repeat-containing / input set: {len(self.repeat_dict)} / {len(self.opmaps)}')
        logger.info(f'Repeat-containing map IDs: {self.repeat_dict.keys()}')
        self.omids_to_draw = list(self.repeat_dict.keys())
        self.repeat_detected = True
        return

    def set_repeat_view(self):
        """Set viewer region by repeat calling results"""
        self.maps_per_page  =   1
        self.reverse        =   False
        self.omids_to_draw  =   sorted(list(self.repeat_dict.keys()))
        curr_map = self.omids_to_draw[self.curr_page-1]
        curr_map_first_repeat_iv = sorted(list(self.repeat_dict[curr_map]))[0]
        curr_map_last_repeat_iv = sorted(list(self.repeat_dict[curr_map]))[-1]
        curr_map_labels     =   self.opmaps.df.loc[curr_map, "labels"]
        logger.debug(f'curr_map_first_repeat_iv: {curr_map_first_repeat_iv}')
        logger.debug(f'curr_map_last_repeat_iv: {curr_map_first_repeat_iv}')
        repeat_start_coord = curr_map_labels[curr_map_first_repeat_iv[0]]
        repeat_end_coord = curr_map_labels[curr_map_last_repeat_iv[1]]
        logger.info(f'repeat_start_coord: {repeat_start_coord}')
        logger.info(f'repeat_last_coord: {repeat_end_coord}')
        view_start_coord = math.floor(repeat_start_coord/1e5) * 1e5
        view_end_coord = max(math.ceil(repeat_end_coord/1e5) * 1e5, view_start_coord + 1e5)
        self.view_region    =   [view_start_coord, view_end_coord]
        logger.debug(f'set_repeat_view -> self.view_region: {self.view_region}')
        logger.debug(f'{repeat_start_coord}, {repeat_end_coord}')
        return

    def get_target_omids(self):
        """Process filtering input for map drawing"""
        # Process omids
        if isinstance(self.omids_to_draw, str):
            if self.omids_to_draw.strip().lower() == 'current':
                logger.debug('Set of maps to draw remains unchanged.')
            elif self.omids_to_draw.strip().lower() == 'all':
                self.omids_to_draw = self.opmaps.df.index
                logger.debug('Including the entire dataset to draw.')
            else:
                tmp_omids_to_draw = []
                omids_to_process = self.omids_to_draw.split(',')
                for x in omids_to_process:
                    try:
                        if ':' in x:
                            start, end = map(int, x.split(':'))
                            tmp_omids_to_draw += range(start, end+1)
                        elif '-' in x:
                            start, end = map(int, x.split('-'))
                            tmp_omids_to_draw += range(start, end+1)
                        else:
                            tmp_omids_to_draw.append(int(x))
                    except ValueError:
                        logger.warning(f'Error in interpreting {x}. {x} is skipped.')
                self.omids_to_draw = tmp_omids_to_draw
                del tmp_omids_to_draw
        elif isinstance(self.omids_to_draw, int):
            self.omids_to_draw = [self.omids_to_draw]
        # validate omids
        self.omids_to_draw  = list(set(self.omids_to_draw))
        tmp_set             = self.omids_to_draw.copy()
        for omid in tmp_set:
            if omid not in self.opmaps.df.index:
                self.omids_to_draw.remove(omid)
                logger.warning(f'Map ID {omid} is not found and is skipped')
        del tmp_set
        # filter repeat if selected
        if self.repeat_only:
            self.omids_to_draw = set(self.repeat_dict).intersection(self.omids_to_draw)
            if self.sort_by == 'id':
                self.omids_to_draw  =   sorted(self.omids_to_draw, reverse=self.reverse)
            elif self.sort_by == 'length':
                self.omids_to_draw  =   sorted(
                                            self.omids_to_draw,
                                            key=lambda x: self.opmaps.df.loc[x, 'labels'][-1],
                                            reverse=self.reverse,
                                        )
            elif self.sort_by == 'repeat_score':
                logger.debug(sorted(list(self.repeat_dict[list(self.omids_to_draw)[0]])))
                self.omids_to_draw  =   sorted(
                                            self.omids_to_draw, 
                                            key=lambda k: 
                                                max(x[2] for x in sorted(list(self.repeat_dict[k]))),
                                            reverse=True,
                                        )
        # sorting
        sorting_order           =   ['forward', 'reverse'][int(self.reverse)]
        logger.debug(f'Sort maps to draw by {self.sort_by} strategy in {sorting_order} order')
        if self.sort_by.lower() == 'id':
            self.omids_to_draw  =   sorted(self.omids_to_draw, reverse=self.reverse)
        elif self.sort_by.lower() == 'length':
            self.omids_to_draw  =   sorted(self.omids_to_draw,
                                            key=lambda x: self.opmaps.df.loc[x, 'labels'][-1],
                                            reverse=self.reverse)
        else:
            self.sort_by        =   'current'
            self.reverse        =   False
        # paginate
        logger.debug(f'get_target_omids->self.curr_page: {self.curr_page}')
        page_start_num          =   (self.curr_page - 1) * self.maps_per_page
        page_end_num            =   page_start_num + self.maps_per_page
        # draw
        omids_on_page           =   self.omids_to_draw[page_start_num:page_end_num]
        logger.info(f'{len(omids_on_page)} map(s) to draw on page: {", ".join(list(map(str, omids_on_page)))}')
        return omids_on_page
    
    def get_total_page(self):
        logger.debug(f'self.omids_to_draw: {self.omids_to_draw}')
        return math.ceil(len(self.omids_to_draw)/self.maps_per_page)

    def draw_maps(  self, omids=None, show_anno=False):
        """Draw opmap objects in ID or length order"""
        fig_rows        =   sum((show_anno, 1))
        row_heights     =   [1]
        if show_anno:
            row_heights =   [0.2, 0.8]
        fig             =   make_subplots(  rows                =   fig_rows,
                                            shared_xaxes        =   True,
                                            vertical_spacing    =   0.02,
                                            row_heights         =   row_heights
                            )
        if show_anno:
            anno_traces =   []              # placeholder
            fig.add_traces(anno_traces)

        map_traces      =   []
        for i, omid in enumerate(omids):
            all_labels          =   self.opmaps.df.loc[omid, 'labels']
            logger.info(f'Map ID: {omid}, Start & End coordinates: {self.view_region}')
            load_region_buffer = 0
            if self.repeat_only:
                logger.info(f'Repeat interval in Map # {omid}: {self.repeat_dict[omid]}')
                logger.info(f'Repeat score for Map # {omid}: {[x[2] for x in sorted(list(self.repeat_dict[omid]))]}')
                highlight_start_coords = []
                highlight_end_coords = []
                logger.debug(f'sorted(self.repeat_dict[omid]): {sorted(self.repeat_dict[omid])}')
                for hightlight_label_pairs in sorted(self.repeat_dict[omid]):
                    logger.debug(f'hightlight_label_pairs: {hightlight_label_pairs}')
                    highlight_start = all_labels[hightlight_label_pairs[0]-1]
                    highlight_start_coords.append(highlight_start)
                    logger.info(f'repeat start: {highlight_start}')
                    highlight_end = all_labels[hightlight_label_pairs[1]-1]
                    logger.info(f'repeat end: {highlight_end}')
                    highlight_end_coords.append(highlight_end)
                first_highlight_start = sorted(highlight_start_coords)[0]
                last_highlight_end = sorted(highlight_end_coords)[-1]
                if first_highlight_start < int(self.view_region[0] - load_region_buffer):
                    self.view_region[0] = int(first_highlight_start/1e5) * 1e5
                if last_highlight_end > int(self.view_region[1] + load_region_buffer):
                    self.view_region[1] = (int(last_highlight_end/1e5)+1) * 1e5
                logger.info(f'highlight_start_coords: {highlight_start_coords}')
                logger.info(f'highlight_end_coords: {highlight_end_coords}')
            else:
                highlight_start_coords    =   None
                highlight_end_coords      =   None
            logger.debug(f'draw_maps.draw_om.self.view_region: {self.view_region}')
            labels = []
            for label in all_labels:
                if label > int(self.view_region[0] - load_region_buffer) and \
                    label < int(self.view_region[1] + load_region_buffer):
                    labels.append(label)
                if label > self.view_region[1] + load_region_buffer:
                    break
            map_traces.append(
                draw_om(
                    omid=omid,
                    labels=labels,
                    y_offset=i,
                    start_coord=max(0, self.view_region[0]),
                    end_coord=min(self.view_region[1], labels[-1]),
                    highlight_start_coords=highlight_start_coords,
                    highlight_end_coords=highlight_end_coords,
                    opmap_length=max(all_labels),
                )
            )
        map_traces          =   list(itertools.chain.from_iterable(map_traces))
        fig.add_traces(map_traces, rows=show_anno + 1, cols=1)
        if len(omids) > 0: # to add auto update view region
            curr_max_len    =   max([len(self.opmaps.df.loc[x]) for x in self.omids_to_draw])
            if curr_max_len >=  self.max_x:
                fig.update_layout(get_layout(
                xaxis       =   {'range': [self.view_region[0], self.view_region[1]]},
                yaxis       =   {'autorange': True},
                height      =   360
                ))
            else:
                fig.update_layout(get_layout(
                    xaxis   =   {'autorange': True},
                    yaxis   =   {'autorange': True},
                    height  =   360
                ))
            if self.maps_per_page <= 10:
                fig.update_layout(get_layout(yaxis={'range':[0,60]}))
        if show_anno:
            fig.update_traces(connectgaps=False)
            fig.update_layout(get_layout(hovermode='x unified'))
        layout = dcc.Graph(figure=fig)
        return layout

    def save_view_as_svg(self):
        """Save view as SVG vector file"""
        return self

    def goto_page_input(self, value: int=1):
        """Input box to set page number"""
        return dcc.Input(   id='page', type='number',
                            placeholder=self.curr_page,
                            value=str(value), min=1, max=self.total_page, step=1)

    def total_page_input(self, total_page):
        """Input box to set total page number; Use Input to retrieve State"""
        return dcc.Input(   id='total_page', type='number',
                            placeholder=total_page,
                            value=str(total_page),
                            min=1, step=0)

    def goto_page_button(self):
        """Diplay custom page of maps"""
        return html.Button( 'Go', id='goto_page', n_clicks=0)

    def next_page_button(self):
        """Diplay next page of maps"""
        return html.Button( '▶', id='next_page', n_clicks=0)

    def prev_page_button(self):
        """Diplay previous page of maps"""
        return html.Button( '◀', id='prev_page', n_clicks=0)

    def maps_per_page_input(self):
        """Diplay number of maps to display per page"""
        return dcc.Input(   id='maps_per_page', type='number',
                            placeholder=10,
                            value=self.maps_per_page,
                            debounce = True)

    def custom_id_input(self):
        """Input custom ID for filtering"""
        omids = ','.join(map(str, self.omids_to_draw))
        return dcc.Input(   id='custom_omids', type='text',
                            placeholder='e.g. 1,2,3,6:12,13-19 OR all',
                            value=omids)
    
    def custom_id_submit(self):
        """Submit custom ID input for filtering"""
        return  html.Button(id='custom_omids_submit',
                            n_clicks=0,
                            children='Go')
    
    def sort_by_button_group(self):
        """Choose map view sorting method"""
        return  dbc.ButtonGroup(
                    id = "sort_by",
                    children = [
                                dbc.Button( 'Map ID', id='by_id', n_clicks=0, active=(self.sort_by=='id')),
                                dbc.Button( 'Length', id='by_length', n_clicks=0, active=(self.sort_by=='length')),
                                dbc.Button( 'Repeat score', id='by_repeat_score', 
                                            n_clicks=0, active=(self.sort_by=='repeat_score'))
                    ]
                )

    def sort_reverse_button(self):
        return  dbc.Button('Reverse', id='reverse', n_clicks=0, active=self.reverse)

    def repeat_only_button(self):
        """Diplay custom page of maps"""
        return  dbc.Button('Repeat Only', id='repeat_only', n_clicks=0, active=self.repeat_only)

    def view_all_button(self):
        """Diplay all maps in dataset"""
        return  dbc.Button('All Maps', id='view_all', n_clicks=0, active=self.view_all)

    def view_region_min_input(self, value: int=None):
        """Input box to set view region min"""
        return dcc.Input(   id='view_region_min', type='number',
                            placeholder=self.view_region[0],
                            value=str(value), min=-1e7, max=3e8, step=5e4,
                            debounce = True)
    
    def view_region_max_input(self, value: int=None):
        """Input box to set view region max"""
        return dcc.Input(   id='view_region_max', type='number',
                            placeholder=self.view_region[1],
                            value=str(value), min=5e4, max=3e8, step=5e4,
                            debounce = True)


def draw_om(omid=None, labels=None,
            start_coord=None, end_coord=None, start_label=None, end_label=None,
            highlight_start_coords=None, highlight_end_coords=None,
            y_offset=1, opmap_length=None, show_label_density=True,
            colors=DEFAULT_COLORSCHEME, color_scheme={}):
    """Draw optical map (Opmap) object as Plotly traces"""
    if len(color_scheme) == 0:
        color_scheme = get_color_scheme(colors)
    if start_coord is None and end_coord is None and \
       start_label is not None and end_label is not None:
        start   =   labels[start_label]
        end     =   labels[end_label]
    else:
        start   =   start_coord
        end     =   end_coord
        if start is None:
            start = 0
        if end is None:
            end = labels[-1]
        start_label = 0

    map_box     =   go.Scatter( x       =   [start, end, end, start, start],
                                y       =   [y - 3 * y_offset for y in [4,4,2,2,4]],
                                name    =   f'Map {omid}',
                                hovertemplate   =   f'Map ID: {omid}' +
                                                    f'<br />Length: {labels[-1]}', #+
                                                    #'<br />Label density: '*show_label_density,
                                fillcolor       =   color_scheme['fillcolor'],
                                fill            =   'toself',
                                marker          =   dict(opacity=0),
                                line_width      =   0,
                                showlegend      =   False
                    )

    highlight_boxes =   []
    if highlight_start_coords is not None:
        assert len(highlight_start_coords) == len(highlight_end_coords)
        for i in range(len(highlight_start_coords)):
            h_start =   highlight_start_coords[i]
            h_end   =   highlight_end_coords[i]
            highlight_boxes.append(
                go.Scatter( x       =   [h_start, h_end, h_end, h_start, h_start],
                            y       =   [y - 3 * y_offset for y in [4,4,2,2,4]],
                            name    =   f'Repeat {omid}_{h_start}_{h_end}',
                            hovertemplate   =   f'<br />Repeat length: {h_end - h_start + 1}', #+
                                                #f'<br />Repeat score: dummy', # to implement
                            fillcolor       =   color_scheme['highlight_color'],
                            fill            =   'toself',
                            marker          =   dict(opacity=0),
                            line_width      =   0,
                            showlegend      =   False
                )
            )
    
    label_lines    =   [go.Scatter( x=[x,x], y=[y - 3 * y_offset for y in [2,4]],
                                    name            =   f'{omid}_{start_label+i}',
                                    hovertemplate   =   f'Map ID: {omid}' +
                                                        f'<br />Label: {start_label+i}' +
                                                        f'<br />Coordinate: {x}' +
                                                        f'<br />Map length: {opmap_length}',
                                    line_color      =   color_scheme['line_color'],
                                    marker          =   dict(opacity=0),
                                    showlegend      =   False
                        ) for i,x in enumerate(labels)]
    map_traces = [map_box] + highlight_boxes + label_lines
    return map_traces


def draw_annotations(anno: Annotations=None):
    colors = px.colors.qualitative.Plotly
    fig = go.Figure()

    features = list(anno.df['feature'].unique())
    for i, row in anno.df.iterrows():
        start   =   row['start']
        end     =   row['end']
        fig.add_traces(
            go.Scatter(
                x = (start, end, end, start, start),
                y=(6,6,5,5,6),
                mode            =   'markers+lines',
                line            =   dict(color=colors[features.index(row['feature'])]),
                name            =   '',
                hovertemplate   =   f'Coordinates: {end}-{start}' +
                                    f'<br />Length: {end - start + 1}<br />' +
                                    '<br />'.join([f'{x.title()}: {row[x]}' for x in ['strand', 'source', 'feature', 'score']]) +
                                    f'<br />{"Attribute"}: {row["attribute"]}',
                fill            =   'toself',
                marker          =   dict(opacity=0),
                line_width      =   0,
                showlegend      =   False
            )
        )
    
    fig.update_traces(connectgaps=False)
    return fig

def repeat_call_settings():
    min_seg_len_setter  =   QuantitySetter( name        =   'min_segment_length',
                                            description =   'Minimum segment length', 
                                            min_val = 50, max_val = 2500, step = 50,
                                            placeholder = 500,
                                            slider_marks = [100, 300, 500, 1000, 2000, 2500]
                            ).layout
    max_k_setter        =   QuantitySetter( name        =   'max_k',
                                            description =   'Maximum K (Number of segments in a repeat unit)', 
                                            min_val = 2, max_val = 31, step = 1,
                                            placeholder = 13,
                                            slider_marks = [3, 13, 23, 31]
                            ).layout
    return html.Div([min_seg_len_setter, max_k_setter])
