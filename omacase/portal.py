#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220324


"""
OMACASE main portal
By Claire Chung
"""


import logging
from argparse import Namespace
from glob import glob
import json
import sys

from dash import dcc, html, Input, Output, State, ALL, callback_context
import dash_bootstrap_components as dbc
import dash_uploader as du

from omacase.app import app
from omacase.config_editor import ConfigEditor, default_args
from omacase.dashboard import ReportDashboard
from omacase.layouts import nav, error_404, upload_box, get_uploaded_dataset_list
from omacase.layouts import dummy_config_html, dummy_omviewer_html, dummy_upload_box
from omacase.omviewer import DisplayMaps
from omacase.annotation_viewer import DisplayAnnotations
from omacase.opmap import check_file_extension, get_opmaps
from omacase.param_parser import parse_opts_from_cmd, UPLOAD_FOLDER_ROOT
from omacase.utils import LOG_FORMAT

logger = logging.getLogger(__name__)


dummy_config_outputs = ['', None, 'success', '']


class Portal(object):
    """OMACASE Homepage"""
    def __init__(self, args=None, app=app):
        if app is not None:
            self._app = app
        if args is not None and hasattr(args, 'verbose') and args.verbose:
            logging.basicConfig(format=LOG_FORMAT, stream=sys.stderr, level=logging.DEBUG)
        self.datasets           =   []
        self.opmaps             =   None
        self.upload_folder_root =   UPLOAD_FOLDER_ROOT
        if self._app is not None and hasattr(self, "callbacks"):
            self.callbacks(self._app)
            self._app.config.suppress_callback_exceptions=True
            du.configure_upload(self._app, UPLOAD_FOLDER_ROOT)
        self.layout             =   self.get_app_layout()
        self.viewer             =   None
        self.conf_editor        =   None
        self.viewer_loaded      =   False
        self.conf_editor_loaded =   False
    
    def callbacks(self, _app):
        @_app.callback( Output('page_content', 'children'),
                        Input('url', 'pathname'),
                        Input({'role': 'dataset', 'index': ALL}, 'n_clicks'),
                        Input({'role': 'dataset', 'index': ALL}, 'value'),
                        State('page_content', 'children'),
                        prevent_initial_call=True) 
        def display_page(
                pathname,
                dataset_n_clicks,
                dataset_name,
                current_layout,
            ):
            changed_id  =   [p['prop_id'] for p in callback_context.triggered][0]
            logger.debug(f'changed_id: {changed_id}')
            split_path  =   pathname.split('/')
            if 'url' in changed_id:
                self.page_object = None
                logger.debug(f'URL path: {pathname}')
                
                if split_path[1] == '':
                    return [None]
                elif split_path[1] == 'qc-report':
                    if self.opmaps is None:
                        return upload_box()
                    else:
                        args = parse_opts_from_cmd()
                        logger.debug(f'args: {args}')
                        changed_id = None
                        return ReportDashboard(opmaps=self.opmaps, app=self._app).layout
                elif split_path[1] == 'map-viewer':
                    args = parse_opts_from_cmd()
                    if len(split_path) >=3 and split_path[2] == 'repeats':
                        args.repeat_only = True
                    else:
                        args.repeat_only = False
                    if self.opmaps is None:
                        return upload_box()
                    else:
                        selected_dataset = self.datasets[-1]
                        logger.debug(check_file_extension(selected_dataset))
                        if check_file_extension(selected_dataset)[0] in ['cmap', 'bnx']:
                            logger.info(f'Loading dataset: {selected_dataset}')
                            self.opmaps = get_opmaps(filepath=selected_dataset, labels_only=True)
                            self.opmaps.filepath = selected_dataset
                            logger.info(f'Loaded dataset: {selected_dataset}')
                        logger.debug(self.opmaps)
                        if not self.viewer_loaded:
                            viewer = DisplayMaps(args=args, opmaps=self.opmaps, app=self._app)
                            self.viewer = viewer
                            self.viewer_loaded = True
                        else:
                            self.opmaps = self.opmaps
                        return self.viewer.layout
                elif split_path[1] == 'annotation-viewer':
                    anno_view = DisplayAnnotations(args=args, q_opmaps=self.opmaps, app=self._app)
                    return anno_view.layout
                elif split_path[1] == 'upload-data':
                    return upload_box()
                elif split_path[1] == 'select-data':
                    return get_uploaded_dataset_list()
                elif split_path[1] == 'config-editor':
                    if not self.conf_editor_loaded:
                        self.conf_editor = ConfigEditor(_app=self._app)
                        self.conf_editor_loaded = True
                    return self.conf_editor.layout
                else:
                    return error_404()
            elif 'upload_data_box' in changed_id:
                return current_layout
            elif 'role' in changed_id and 'dataset' in changed_id:
                i = json.loads(changed_id.split('.')[0])['index']
                selected_dataset = sorted(glob(f'{self.upload_folder_root}*/*'))[i]
                logger.info(f'Selected dataset: {selected_dataset}')
                self.datasets.append(selected_dataset)
                selected_dataset = self.datasets[-1]
                logger.debug(check_file_extension(selected_dataset))
                if check_file_extension(selected_dataset)[0] in ['cmap', 'bnx']:
                    logger.info(f'Loading dataset: {selected_dataset}')
                    self.opmaps = get_opmaps(filepath=selected_dataset, stat_only=True)
                    self.opmaps.filepath = selected_dataset
                    logger.info(f'Loaded dataset: {selected_dataset}')
                return current_layout

    def get_app_layout(self):
        return  html.Div(
                    [
                        nav(),
                        dcc.Location(id='url', href='/', refresh=False),
                        dcc.Store(id='memory_output', storage_type='local'),
                        dcc.Loading(    
                            id='page_content',
                            className='container',
                            children=[
                                html.Div(
                                    className='row',
                                    children=html.H3(
                                    className='col-12',
                                    children='Please select a menu item to start',
                                    ),
                                ),
                            ],
                        ),
                        dummy_config_html,
                        dummy_omviewer_html,
                        dummy_upload_box,
                    ],
                )


    def start(self):
        self._app.layout =  html.Div([
                                self._app.layout.children[0],
                                self.layout,
                                self._app.layout.children[1]
                            ])
        self._app.run_server(host='0.0.0.0', port=args.PORT, debug=True, use_reloader=False)

if __name__ == "__main__":
    args = Namespace(**default_args)
    Portal(args=args, app=app).start()
