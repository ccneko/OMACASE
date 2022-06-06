#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220516


"""
OMACASE, a tool to analyze optical mapping data
with _de novo_ repeat motif detection and analysis capabilities
distributed under license GPL-3.0-or-later
By Claire Chung
"""


import logging
from sys import argv, stderr
from unittest.mock import MagicMock

from dash import html
from gevent.pywsgi import WSGIServer
import pdfkit

from omacase import __path__ as project_paths
from omacase import utils, param_parser, report, repeat_caller, repeat_simulator
from omacase.utils import LOG_FORMAT
from omacase import dashboard, omviewer, portal
from omacase.app import app


logging.basicConfig(format=LOG_FORMAT, stream=stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def annotation_viewer_mode(args):
    raise NotImplementedError



def PDFKit_patched(
        url_or_file, type_, options, toc, cover,
        css, configuration, cover_first, verbose, 
        mock_handle_error,
    ):
    """Patch PDFKit error handling to bypass wkhtmltopdf ContentNotFound error"""
    mock_handle_error.return_value = None
    #return  pdfkit.PDFKit(
                #url_or_file=url_or_file, type_=type_, options=options, toc=toc, 
                #cover=cover, css=css, configuration=configuration,
                #cover_first=cover_first, verbose=verbose
    #        )



FUNCTION_MAP    =   {
                        'text-report':  report.text_report,
                        'report':       dashboard.ReportDashboard,
                        'repeat-bed':   repeat_caller.write_repeat_bed,
                        'repeat-sim':   repeat_simulator.repeat_sim,
                        'view':         omviewer.DisplayMaps,
                        'web' :         portal.Portal,
                        'anno':         annotation_viewer_mode,\
                        'p':            report.text_report,
                        'r':            dashboard.ReportDashboard,
                        'b':            repeat_caller.write_repeat_bed,
                        'v':            omviewer.DisplayMaps,
                        'w' :           portal.Portal,
                        'a':            annotation_viewer_mode,
                        's':            repeat_simulator.repeat_sim,
                    }



if __name__ == "__main__":
    utils.banner_art()
    logger.debug(argv)
    args = param_parser.parse_opts_from_cmd(argv)
    logger.info(vars(args))
    
    if args.mode.lower() in FUNCTION_MAP:
        func = FUNCTION_MAP[args.mode.lower()]
        if args.mode.lower() in ['text-report', 'repeat-bed', 'repeat-sim']:
            func(args)
        else:
            logger.info(vars(args))
            if not hasattr (args, 'report_output_format') \
                or args.report_output_format is None \
                or 'pdf' not in args.report_output_format:
                app.layout  =   html.Div([
                                app.layout.children[0],
                                func(args=args, app=app).layout,
                                app.layout.children[1]
                            ])
                logger.info('Server started')
                WSGIServer(('0.0.0.0', args.port), app.server).serve_forever()
            elif func == dashboard.ReportDashboard:
                html_out = func(args).layout
                css_folder = project_paths[0] + '/assets/css/'
                stylesheets = [f'{css_folder}00_bootstrap.min.css', f'{css_folder}01_omacase.css']
                if args.output_prefix is not None:
                    output_prefix = args.output_prefix
                else:
                    output_prefix = f'reports/{utils.timestamp()}-omacase-output'
                output_path=f'{output_prefix}.pdf'
                PDFKit_patched  =   pdfkit.PDFKit(
                                        url_or_file=str(html_out),
                                        type_='string',
                                        css=stylesheets
                                    )
                PDFKit_patched.handle_error = MagicMock(return_value=None, stderr_lines='')
                PDFKit_patched.to_pdf(output_path)