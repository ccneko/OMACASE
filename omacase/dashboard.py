#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220331


"""
Dashboard OM dataset report webpage components
By Claire Chung
"""


import logging

from dash import html
from holoviews.plotting.plotly.dash import to_dash

from omacase.opmap import Opmaps, get_opmaps
from omacase.report import opmaps_summary_dict, opmaps_summary_df
from omacase.dataviz import *
from omacase.layouts import dash_df
from omacase.layouts import dummy_config_html, dummy_omviewer_html
from omacase.utils import render_template_without_request



logger = logging.getLogger(__name__)


class ReportDashboard(object):
    """
    Dashboard web app for Opmaps OM dataset report
    """
    def __init__(self, args=None, opmaps: Opmaps=None, app=None, interactive=False):
        """Initializes ReportDashboard object"""
        self.interactive = interactive
        if opmaps is not None:
            if not isinstance(opmaps, Opmaps):
                raise ValueError(f'Input (type {type(opmaps)}) is not Opmaps OM dataset object.')
            self.opmaps         =   opmaps
        else:
            self.opmaps         =   get_opmaps(args=args, stat_only=True)
        if not hasattr (args, 'report_output_format') \
            or args.report_output_format is None \
            or 'pdf' not in args.report_output_format:
            self.layout         =   self.get_app_layout(app)
        else:
            self.layout     =   self.get_static_html()
    
    def get_app_layout(self, app=None):
        """Set OMACASE dashboard web app layout"""
        if self.opmaps.file_format == 'bnx':
            qx_plots_results = qx_plots(self.opmaps, interactive=self.interactive)
        if self.interactive:
            if self.opmaps.file_format == 'bnx':
                qx11_plot_html      =   html.Div(
                                            id = 'qx11_hist',
                                            children = to_dash(
                                                app,
                                                qx_plots_results[0]
                                        ).children)
                qx12_plot_html      =   html.Div(
                                            id = 'qx12_hist',
                                            children = to_dash(
                                                app,
                                                qx_plots_results[1]
                                        ).children)
                mol_snr_plot_html   =   html.Div(
                                            id = 'mol_snr_hist',
                                            children = to_dash(
                                                app,
                                                qx_plots_results[2]
                                        ).children)
                mol_inten_plot_html   =   html.Div(
                                            id = 'mol_intensity_hist',
                                            children = to_dash(
                                                app,
                                                qx_plots_results[3]
                                        ).children)
            map_len_hist_html   =   html.Div(
                                        id = 'map_len_hist',
                                        children = to_dash(
                                            app,
                                            map_len_hist(
                                                self.opmaps,
                                                interactive=self.interactive,
                                            )
                                    ).children)
            n_metrics_hist_html =   html.Div(
                                        id = 'n_metrics_hist',
                                        children = to_dash(
                                            app,
                                            n_metrics_hist(self.opmaps)
                                    ).children)
        else:
            if self.opmaps.file_format == 'bnx':
                qx11_plot_html      =   html.Img(
                                            id = 'qx11_plot',
                                            src=qx_plots_results[0],
                                        )
                qx12_plot_html      =   html.Img(
                                            id = 'qx12_plot',
                                            src=qx_plots_results[1],
                                        )
                mol_snr_plot_html   =   html.Img(
                                            id = 'mol_snr_plot',
                                            src=qx_plots_results[2],
                                        )
                mol_inten_plot_html =   html.Img(
                                            id = 'mol_intensity_plot',
                                            src=qx_plots_results[3],
                                        )
            map_len_hist_html       =   html.Img(
                                            id = 'map_len_hist',
                                            src=map_len_hist(
                                                self.opmaps,
                                                interactive=self.interactive,
                                            ),
                                        )
            n_metrics_hist_html     =   html.Img(
                                            id = 'n_metrics_hist',
                                            src=n_metrics_hist(
                                                self.opmaps,
                                                interactive=self.interactive,
                                            ),
                                        )
        # body sections
        report_table =  html.Div(
                            className='row',
                            children=[
                                html.Div(
                                    className='col-12',
                                    children=[
                                        html.H2(children='Dataset summary'),
                                        html.Div(id='opmaps_summary_div', children=[
                                                dash_df(opmaps_summary_df(self.opmaps),
                                                    'opmaps_summary')]),
                                    ],
                                ),
                            ],
                        )

        if self.opmaps.file_format == 'bnx':
            plots_section = html.Div(id='plots', children=[
                                html.Div(
                                    className='row',
                                    children=[
                                        html.H2(
                                            className='col-12 col-md-4',
                                            children='Map length statistics',
                                        ),
                                        html.H2(
                                            className='col-12 col-md-4',
                                            children='Map quality statistics (Labels)',
                                        ),
                                        html.H2(
                                            className='col-12 col-md-4',
                                            children='Map quality statistics (Mol.)',
                                        ),
                                    ]
                                ),
                                html.Div(
                                    className='row',
                                    children=[
                                        html.H3(
                                            className='col-12 col-md-4',
                                            children='Map length histogram',
                                        ),
                                        html.H3(
                                            className='col-12 col-md-4',
                                            children='Per molecule label SNR',
                                        ),
                                        html.H3(
                                            className='col-12 col-md-4',
                                            children='Per molecule backbone SNR'
                                        ),
                                    ]
                                ),
                                html.Div(
                                    className='row',
                                    children=[
                                        html.P(
                                            className='col-12 col-md-4',
                                            children='Pink: 20-150 kbp; Yellow: 150-220 kbp; Green: >220 kbp',
                                        ),
                                        html.Div(className='col-12 col-md-8'),
                                    ],
                                ),
                                html.Div(
                                    className='row',
                                    children=[
                                        html.Div(
                                            className='col-12 col-md-4',
                                            children=map_len_hist_html,
                                        ),
                                        html.Div(
                                            className='col-12 col-md-4',
                                            children=qx11_plot_html,
                                        ),
                                        html.Div(
                                            className='col-12 col-md-4',
                                            children=mol_snr_plot_html,
                                        ),
                                ]),
                                html.Div(
                                    className='row',
                                    children=[
                                        html.H3(
                                            className='col-12 col-md-4',
                                            children='Map length N-metrics',
                                        ),
                                        html.H3(
                                            className='col-12 col-md-4',
                                            children='Per molecule label intensity',
                                        ),
                                        html.H3(
                                            className='col-12 col-md-4',
                                            children='Per molecule backbone intensity',
                                        ),
                                    ],
                                ),
                                html.Div(
                                    className='row',
                                    children=[
                                    html.Div(
                                        className='col-12 col-md-4',
                                        children=n_metrics_hist_html,
                                    ),
                                    html.Div(
                                        className='col-12 col-md-4',
                                        children=qx12_plot_html,
                                    ),
                                    html.Div(
                                        className='col-12 col-md-4',
                                        children=mol_inten_plot_html,
                                    ),
                                ]),
                            ])
        else:
            plots_section = html.Div(
                                id='plots',
                                className='row',
                                children=[
                                    html.Div(
                                        className='col-12',
                                        children=[
                                            html.H2(children='Map length statistics'),
                                            html.H3(children='Map length histogram'),
                                            html.P(children='Pink: 20-150 kbp; Yellow: 150-220 kbp; Green: >220 kbp'),
                                            map_len_hist_html,
                                            html.H3(children='Map length N-metrics'),
                                            n_metrics_hist_html,
                                        ],
                                    ),
                                ],
                            )

        return  html.Div(
                    id='report_body',
                    className='row page-content',
                    children=[
                        html.Div(className='col-1'),
                        html.Div(
                            className='col-10',
                            children=[
                                report_table,
                                html.Div(className='row sep-row'),
                                plots_section,
                                dummy_config_html,
                        ]),
                        html.Div(className='col-1'),
                    ],
                )


    def get_static_html(self):
        """Set OMACASE dashboard static HTML layout"""
        qx_plots_results    = qx_plots(self.opmaps, width=360, height=340)
        map_len_hist_src    = map_len_hist(self.opmaps, width=360, height=340)
        n_metrics_hist_src  = n_metrics_hist(self.opmaps, width=360, height=340)
        om_info_dict        = opmaps_summary_dict(self.opmaps)
        if self.opmaps.file_format == 'bnx':
            template = 'report-bnx.html'
            layout = render_template_without_request(
                template,
                time_now            =   om_info_dict['Report generation datetime'],
                filename            =   om_info_dict['Filename'],
                file_format         =   om_info_dict['File type'],
                file_version        =   om_info_dict['File format version'],
                motif               =   om_info_dict['Label motif'],
                map_num             =   om_info_dict['Map number'],
                total_length        =   om_info_dict['Total length'],
                n50                 =   om_info_dict['N50 length (bp)'],
                l50                 =   om_info_dict['L50 value'],
                n90                 =   om_info_dict['N90 length (bp)'],
                l90                 =   om_info_dict['L90 value'],
                max_length          =   om_info_dict['Max length (bp)'],
                min_length          =   om_info_dict['Min length (bp)'],
                label_density       =   om_info_dict['Label density / 100 kbp'],

                qx11_plot_src       =   qx_plots_results[0],
                qx12_plot_src       =   qx_plots_results[1],
                mol_snr_plot_src    =   qx_plots_results[2],
                mol_inten_plot_src  =   qx_plots_results[3],
                map_len_hist_src    =   map_len_hist_src,
                n_metrics_hist_src  =   n_metrics_hist_src,
                qx11_avg            =   om_info_dict['Mean of mol-average label SNR'],
                qx11_sd             =   om_info_dict['SD of mol-average label SNR'],
                qx12_avg            =   om_info_dict['Mean of mol-average label intensity'],
                qx12_sd             =   om_info_dict['SD of mol-average label intensity'],
                mol_snr_avg         =   om_info_dict['Mean of mol-average mol SNR'],
                mol_snr_sd          =   om_info_dict['SD of mol-average mol SNR'],
                mol_intensity_avg   =   om_info_dict['Mean of mol-average mol intensity'],
                mol_intensity_sd    =   om_info_dict['SD of mol-average mol intensity'],
            )
        elif self.opmaps.file_format == 'cmap':
            template = 'report-cmap.html'
            layout = render_template_without_request(
                    template,
                    time_now            =   om_info_dict['Report generation datetime'],
                    filename            =   om_info_dict['Filename'],
                    file_format         =   om_info_dict['File type'],
                    file_version        =   om_info_dict['File format version'],
                    motif               =   om_info_dict['Label motif'],
                    map_num             =   om_info_dict['Map number'],
                    total_length        =   om_info_dict['Total length'],
                    n50                 =   om_info_dict['N50 length (bp)'],
                    l50                 =   om_info_dict['L50 value'],
                    n90                 =   om_info_dict['N90 length (bp)'],
                    l90                 =   om_info_dict['L90 value'],
                    max_length          =   om_info_dict['Max length (bp)'],
                    min_length          =   om_info_dict['Min length (bp)'],
                    label_density       =   om_info_dict['Label density / 100 kbp'],
                    map_len_hist_src    =   map_len_hist_src,
                    n_metrics_hist_src  =   n_metrics_hist_src,
                )
        logger.info('Finish layout rendering')
        return layout
