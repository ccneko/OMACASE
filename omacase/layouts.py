#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220324


"""
Webpage UI components
By Claire Chung
"""

import logging
from datetime import datetime
from glob import glob

from dash import dcc, html, Input, Output, dash_table
import dash_bootstrap_components as dbc
import dash_uploader as du

from omacase.param_parser import MAX_DATASETS, UPLOAD_FOLDER_ROOT

logger = logging.getLogger(__name__)


def web_header():
    return  html.Header(
                id='omacase-header',
                children=[   
                    html.Div([
                        html.H1 (id='home_logo', children=html.A('OMACASE', id='home_logo_link', href='/')),
                        html.P  (id='slogan', children="Optical mapping QC & repeat calling toolkit")
                    ]),
                ],
    )

def web_footer():
    """ReportDashboard footer"""
    copyright_statement =   'Copyright © ' + str(datetime.now().year) + ' Claire Chung'
    cuhk_link           =   html.A( 'CUHK', id='cuhk_link', href='https://www.cuhk.edu.hk')
    ominfo_link         =   html.A( 'OMInfo', id='ominfo_link',
                                    href='https://www.opticalmapping.info')
    cite_link           =   html.Cite(html.A( 'Cite', id='cite_link',
                                    href='/dummy-cite-link'))
    footer              =   html.Footer(
                                id='omacase-footer',
                                children = [
                                    copyright_statement,
                                    cuhk_link,
                                    ominfo_link,
                                    cite_link,
                                ],
                            )
    return  footer

def nav():
    return  dbc.Nav([
                dbc.NavItem(dbc.NavLink(id='upload_data', 
                    children='Upload Dataset', href='/upload-data')),
                dbc.NavItem(dbc.NavLink(id='select_data', 
                    children='Select Uploaded Dataset', href='/select-data')),
                dbc.NavItem(dbc.NavLink( id='qc_report_link',
                    children='Dataset QC', href='/qc-report')),
                dbc.NavItem(dbc.NavLink( id='display_maps_link',
                    children='Map Viewer', href='/map-viewer')),
                #dbc.NavItem(dbc.NavLink( id='display_annotations_link',
                #    children='Display Annotation', href='/annotation-viewer')),
                dbc.NavItem(dbc.NavLink( id='edit_config_link',
                    children='Edit Config', href='/config-editor')),
            ])


def get_uploaded_dataset_list(folder=UPLOAD_FOLDER_ROOT):
        return  html.Div(
                    id='uploaded-list-box',
                    className='row',
                    children=[
                        html.Div(className='col-1'),
                        html.Div(
                            className='col-10',
                            children=[
                                html.H2('Uploaded Datasets'),
                                html.Ul(
                                    id='uploaded-list',
                                    className='uploaded-list',
                                    children=[
                                        html.Li(
                                            dbc.Button(
                                                id={'role': 'dataset', 'index': i},
                                                active=False,
                                                children=x
                                        ),
                                        ) for i,x in enumerate(sorted(glob(f'{folder}*/*')))
                                    ]
                                ),
                            ],
                        ),
                        html.Div(className='col-1'),
                        dummy_config_html,
                    ],
                )

class QuantitySetter:
    def __init__(   self, name, description, min_val=0, max_val=None, step=1,
                    placeholder=0, slider_marks=[]
                ):
        self.name               =   name
        self.description        =   description
        self.min                =   min_val
        self.max                =   max_val
        self.step               =   step
        self.placeholder        =   placeholder
        self.marks              =   slider_marks
        self.layout             =   html.Div([  html.H3(self.description),
                                                self.input(self.placeholder),
                                                self.slider(self.placeholder)])

    def quantity_setter_callbacks(self, app=None):
        app.callback(
                        Output(self.name + '_quantity_input', 'value'),
                        Input(self.name + '_slider', 'value'),
                        prevent_initial_call=True
        )(str)
        app.callback(
                        Output(self.name + '_slider', 'value'),
                        Input(self.name + '_input', 'value'),
                        prevent_initial_call=True
        )(int)
        return  self

    def input(self, value: int=None):
        """Input box to set repeat calling resolution (bp)"""
        return  dcc.Input(      id          =   self.name + '_input',
                                type        =   'number', 
                                placeholder =   self.placeholder,
                                value       =   str(value))

    def slider(self, value: int=None):
        """Slider to set repeat calling resolution (bp)"""
        return  dcc.Slider(
                                id      =   self.name + '_slider',
                                min     =   self.min,
                                max     =   self.max,
                                step    =   self.step,
                                marks   =   dict(zip(   list(map(str,self.marks)),
                                                        self.marks)),
                                value   =   int(value),
                    )


def choose_output_config_format():
    """Choose config output format menu"""
    return  [   
                html.H3(        children='Output config file format'),
                dcc.Checklist(  id='config_format_choice',
                                options=[
                                    {'label': 'JSON', 'value': 'json'},
                                    {'label': 'TOML', 'value': 'toml'}
                                ],
                                value=['json']
                ),
            ]


def slide_out_menu(components=[]):
    """Slide out OMACASE config menu"""
    components      =   []
    menu_components =   html.Div(id='slide', children=[html.Div(id='toggle')] + components)
    logger.debug('menu: {}'.format(menu_components))
    return  menu_components


def dash_df(df, df_id):
    """Create dash_table.DataTable from pandas dataframe"""
    return  dash_table.DataTable(
                    id      =   df_id,
                    columns =   [{'name': i, 'id': i} for i in df.columns],
                    data    =   df.to_dict('records'),
            )


def upload_box():
    return  html.Div(
                className='container',
                children=[
                    du.Upload(  
                        id='upload_data_box',
                        filetypes=['gz', 'bnx', 'cmap', 'gtf'],
                        text='Drag and Drop files or Click here',
                        text_completed='Upload complete: ',
                        pause_button=False,
                        cancel_button=True,
                        chunk_size = 100, 
                        max_file_size=10240,
                    ),
                    dummy_config_html,
                    dummy_omviewer_html,
                ],
            )


def error_404():
    """Error 404 layout"""
    return  html.Div(    id='404',
                        children=[
                            html.H1('404 Not found'),
                            html.P('Please click on one of the links in the navigation menu.')]
            )


dummy_config_html = html.Div(
    className='hidden',
    children=[
        html.Div(id='dummy1', className='hidden'),
        html.Div(id='dummy2', className='hidden'),
        html.Div(id='dummy3', className='hidden'),
        html.Div(id='dummy4', className='hidden'),
        html.Div(id='export_config_btn', className='hidden'),
        html.Div(id='import_config_btn', className='hidden'),
        dcc.Textarea(id='config_format_selector', className='hidden', value=''),
        dcc.Textarea(id='config_editor', className='hidden', value=''),
        dbc.Alert(
            id='config_validator', className='hidden',
            children='', color='success',
        ),
        dcc.Download(id='export_config', data=None),
    ],
)

dummy_omviewer_html = html.Div(
    className='hidden',
    children=[
        dcc.Input(id='page', value=''),
        dcc.Input(id='total_page', value=''),
        dcc.Input(id='custom_omids', value=''),
        html.Button(id='custom_omids_submit', n_clicks=0, value=''),
        dbc.Button(id='by_id', n_clicks=0, active=False),
        dbc.Button(id='by_length', n_clicks=0, active=False),
        dbc.Button(id='by_repeat_score', n_clicks=0, active=False),
        dbc.Button(id='reverse', n_clicks=0, active=False),
        dbc.Button(id='repeat_only', n_clicks=0, active=False),
        dbc.Button(id='view_all', n_clicks=0, active=False),
        dcc.Input(id='view_region_min', type='number', value=0, debounce=True),
        dcc.Input(id='view_region_max', type='number', value=0, debounce=True),
        html.Div(id='omviewer'),
        html.Button(id='prev_page', n_clicks=0),
        html.Button(id='next_page', n_clicks=0),
        html.Button(id='goto_page', n_clicks=0),
        dcc.Input(id='maps_per_page', type='number', value=0, debounce=True),
    ]
)

dummy_upload_box = html.Div(id='upload_data_box')