import logging
import glob

from dash import Dash, html

from omacase import __path__ as project_paths
from omacase.layouts import web_header, web_footer

logger = logging.getLogger(__name__)

external_stylesheets = sorted(glob.glob(f'{project_paths[0]}/assets/css/*.css', recursive=True))
external_scripts = sorted(glob.glob(f'{project_paths[0]}/assets/css/*.css', recursive=True))
logger.debug(f'Stylesheets loaded: {external_stylesheets}')

app         =   Dash(__name__,
                    external_stylesheets=external_stylesheets,
                    external_scripts=external_scripts
                )
app.layout  =   html.Div([web_header(), web_footer()])
app.title   =   'OMACASE: Optical Map Toolkit'
