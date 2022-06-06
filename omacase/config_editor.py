#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220324


"""
OMACASE main portal
By Claire Chung
"""


import logging
import json

import toml
from dash import dcc, html, Input, Output, State, callback_context
import dash_bootstrap_components as dbc

from omacase import __path__ as project_paths
from omacase.app import app
from omacase.layouts import dummy_omviewer_html

logger = logging.getLogger(__name__)

default_config_path =   f'{project_paths[0]}/config/default.toml'
default_args        =   toml.load(default_config_path)

class ConfigEditor(object):
    def __init__(self, format='toml', _app=app):
        if _app is not None:
            self._app = app
        config = ''.join(open(default_config_path).readlines())
        if format == 'json':
            config = json.dumps(config, indent=4)
        format_selector =   dcc.Dropdown(   id          =   'config_format_selector',
                                            options     =   [
                                                                {'label': 'json', 'value': 'json'},
                                                                {'label': 'toml', 'value': 'toml'}
                                                            ],
                                            value       =   'toml')
        editor          =   dcc.Textarea(   id          =   'config_editor',
                                            value       =   config,
                                            autoFocus   =   'autoFocus',
                                            cols        =   100,
                                            className   =   '.mt-2 .mb-2')
        
        validator       =   dbc.Alert(      id          =   'config_validator', children='Valid format', 
                                            color       =   'success',
                                            dismissable =   False,
                                            fade        =   False)
        export_button   =   html.Div([      html.Button(
                                                "Validate & Export",
                                                id="export_config_btn",
                                                className="col-2",
                                            ),
                                            dcc.Download(id="export_config")])
        import_button   =   html.Div([      html.Button(
                                                "Validate & Import",
                                                id="import_config_btn",
                                                className="col-2",
                                            )])    
        self.layout     =   html.Div(
                                id='config_content',
                                className='page-content',
                                children=[
                                    html.Div(
                                        className='row',
                                        children=[
                                            html.Div(className='col-1'),
                                            html.Div(
                                                className='col-10',
                                                children=[
                                                    format_selector, editor, validator,
                                                ],
                                            ),
                                            html.Div(className='col-1'),
                                        ],
                                    ),
                                    html.Div(
                                        className='row',
                                        children=[
                                            html.Div(className='col-1'),
                                            export_button,
                                            import_button,
                                            html.Div(className='col-1'),
                                        ]
                                    ),
                                    dummy_omviewer_html,
                                ],
                            )
        self.callbacks(self._app)

    def callbacks(self, _app):
        @_app.callback( Output('config_editor', 'value'),
                        Output('config_validator', 'children'),
                        Output('config_validator', 'color'),
                        Output('export_config', 'data'),
                        Input('export_config_btn', 'n_clicks'),
                        Input('import_config_btn', 'n_clicks'),
                        Input('config_format_selector', 'value'),
                        State('config_editor', 'value'),
                        State('config_format_selector', 'value'),
                        prevent_initial_call=True)
        def update_config_validator(
                config_export_btn_n, config_import_btn_n, 
                to_format, config, config_format
            ):
            changed_id  =   [p['prop_id'] for p in callback_context.triggered][0]
            logger.debug(f'Last changed element ID: {changed_id}')
            if 'export_config' in changed_id or 'import_config' in changed_id:
                try:
                    if config_format == 'json':
                        args    =   json.loads(config)
                        config  =   json.dumps(args, indent=4)
                    elif config_format == 'toml':
                        args    =   toml.loads(config)
                        config  =   toml.dumps(args)
                    invalid_params  =   set(args) - set(default_args)
                    if len(invalid_params) == 0:
                        return  config, 'Valid format', 'success',\
                                dict(content=config, filename=f'config.{config_format}')
                    else:
                        config_status  =   f'Invalid parameter(s): {invalid_params}'
                        return config, config_status, 'danger', None
                except Exception as e:
                    config_status      =   f'Invalid format. Error message: {e}'
                    return config, config_status, 'danger', None
            elif changed_id == 'config_format_selector.value':
                try: 
                    if to_format == 'json':
                        args        =   toml.loads(config)
                        config      =   json.dumps(args, indent=4)
                    elif to_format == 'toml':
                        args        =   json.loads(config)
                        config      =   toml.dumps(args)
                    return config, f'Changed to {to_format.upper()} format', 'info', None
                except Exception as e:
                    config_status   =   f'Invalid operation. Error message: {e}'
                    return config, config_status, 'danger', None


if __name__ == "__main__":
    ConfigEditor()