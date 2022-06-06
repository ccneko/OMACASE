#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220219


"""
Generate OM summary data visualizations
By Claire Chung
"""


import logging
import io
import base64

import numpy as np
import holoviews as hv
from holoviews import opts
from holoviews.plotting.plotly.dash import to_dash
import plotly.graph_objects as go
import plotly.express as px
from  plotly.subplots import make_subplots
from omacase.opmap import Opmaps, BNX, RMAP
from omacase.alignment import OMAlignments


logger = logging.getLogger(__name__)


def n_metrics_hist(opmaps: Opmaps, interactive=False, width=400, height=400):
    """Plot N10 ~ N90 histogram for Opmaps object"""
    plot = hv.Bars([(f'N{i}', opmaps.n_metrics(i)[0]) for i in list(range(10,101,10))],\
                    'N-metrics', 'Length (bp)')
    if interactive:
        return  plot
    else:
        buf = io.BytesIO()
        if opmaps.file_format == 'cmap':
            width = 1200
            height = 400
        plot = plot.opts(width=width, height=height)
        hv.save(plot, buf, fmt='png')
        data = base64.b64encode(buf.getbuffer()).decode("utf8")
        return "data:image/png;base64,{}".format(data)


def map_len_hist(
        opmaps: Opmaps, interactive=False, bins=500, ideal_n50=220000,
        width=400, height=400, 
    ):
    """Opmaps map length histogram panel"""
    if opmaps is None:
        logger.critical('No Opmaps object loaded.')
        return
    n50             =   opmaps.get_n50()[0]
    plot            =   hv.Histogram(np.histogram(opmaps.lengths, bins=bins), vdims='Map number')\
                        .opts(color='navy', line_color='navy')
    n50_color       =   ['#FF0080', '#27E665'][n50>=ideal_n50]
    n50_vline       =   hv.VLine(n50).opts(opts.VLine(line_color=n50_color, line_width=3,
                        labelled=['N50', n50]))
    #n50_label       =   hv.Text(n50, 20, f'N50: {n50:,}')
    qc_220k_vline   =   hv.VLine(ideal_n50).opts(opts.VLine(line_color='#29DFF0', line_width=3,
                        labelled=['Ideal min N50', ideal_n50]))
    if opmaps.file_format == 'bnx':
        no_vspan    =   hv.VSpan(0, 20000).opts(fillcolor='gray', line_color='gray')
        ng_vspan    =   hv.VSpan(20000, 150000).opts(fillcolor='tomato', line_color='tomato')
        ok_vspan    =   hv.VSpan(150000, ideal_n50).opts(fillcolor='yellow', line_color='yellow')
        gd_vspan    =   hv.VSpan(ideal_n50, 3000000).opts(fillcolor='green', line_color='green')
        plot        =   no_vspan * ng_vspan * ok_vspan * gd_vspan * plot * qc_220k_vline \
                        * n50_vline.opts(xlabel='Length (bp)')
    else:
        plot        =   plot * n50_vline
    if interactive:
        return plot
    else:
        buf = io.BytesIO()
        if opmaps.file_format == 'cmap':
            width   = 1200
            height  = 400    
        plot = plot.opts(width=width, height=height)
        hv.save(plot, buf, fmt='png')
        data = base64.b64encode(buf.getbuffer()).decode("utf8")
        return f'data:image/png;base64,{data}'


def repeat_arr_len_hist(rmap_filepath=None, width=400, height=400):
    """Repeat array length histogram panel"""
    if rmap_filepath == None:
        rmap_filepath = '../data/SAMN02143031-S_Up_H32_molecules.bnx-rep-solve-3.6.out.rmap'
    df          =   RMAP(rmap_filepath).df
    dataset     =   hv.Dataset(df)
    plot        =   hv.operation.histogram(dataset, dimension='AvgRepeatLength', normed=False)\
                    .opts( 
                        hooks=[pan_drag], bgcolor='white',
                        line_width=1, line_color="aliceblue",
                        width=width, height=height,
                    )
    return plot


def repeat_cn_hist(rmap_filepath=None, width=400, height=400):
    """Repeat array copy number histogram panel"""
    if rmap_filepath == None:
        rmap_filepath = '../data/SAMN02143031-S_Up_H32_molecules.bnx-rep-solve-3.6.out.rmap'
    df          =   RMAP(rmap_filepath).df
    dataset     =   hv.Dataset(df)
    plot        =   hv.operation.histogram(dataset, dimension='NRepeatUnits', normed=False)\
                    .opts( 
                        hooks=[pan_drag], bgcolor='white', 
                        line_width=1, line_color="aliceblue",
                        width=width, height=height,
                    )
    return plot

        
def qx_plots(
        opmaps: BNX,
        draw_for=['qx11', 'qx12', 'mol_snr', 'mol_intensity'],
        interactive=False,
        width=400, height=400,
        ):
    if opmaps.file_format != 'bnx':
        return None
    
    plots = []
    if  'qx11'          in draw_for:
        metrics     =   'qx11_mol_mean'
        per_mol_qx  =   opmaps.df[metrics]
        ylabel    =   'Per molecule label SNR'
        plot_color  =   'lightseagreen'
        plots.append(hv.Violin(per_mol_qx, vdims=metrics)\
            .opts(opts.Violin(color=plot_color, width=width, height=height, ylabel=ylabel)))
    if  'qx12'          in draw_for:
        metrics     =   'qx12_mol_mean'
        per_mol_qx  =   opmaps.df[metrics]
        ylabel    =   'Per molecule label intensity'
        plot_color  =   'darkgreen'
        plots.append(hv.Violin(per_mol_qx, vdims=metrics)\
            .opts(opts.Violin(color=plot_color, width=width, height=height, ylabel=ylabel)))
    if  'mol_snr'       in draw_for:
        metrics     =   'mol_snr'
        per_mol_qx  =   opmaps.df[metrics]
        ylabel    =   'Per molecule molecule SNR'
        plot_color  =   'skyblue'
        plots.append(hv.Violin(per_mol_qx, vdims=metrics)\
            .opts(opts.Violin(color=plot_color, width=width, height=height, ylabel=ylabel)))
    if  'mol_intensity' in draw_for:
        metrics     =   'mol_intensity'
        per_mol_qx  =   opmaps.df[metrics]
        ylabel    =   'Per molecule molecule intensity'
        plot_color  =   'slateblue'
        plots.append(hv.Violin(per_mol_qx, vdims=metrics)\
            .opts(opts.Violin(color=plot_color, width=width, height=height, ylabel=ylabel)))
    if interactive:
        return plots
    else:
        pngs = []
        for plot in plots:
            buf = io.BytesIO()
            hv.save(plot, buf, fmt='png')
            data = base64.b64encode(buf.getbuffer()).decode("utf8")
            pngs.append(f'data:image/png;base64,{data}')
        return pngs


def draw_anno(annotations=None, format='gtf', alignments=None, r_opmaps=None, q_opmaps=None):
    colors = px.colors.qualitative.Plotly
    fig = make_subplots(rows=2, shared_xaxes=True, vertical_spacing=0.02)

    features = list(annotations.df['feature'].unique())
    if annotations is not None:
        for i, row in annotations.df.iterrows():
            start   =   row['start']
            end     =   row['end']
            fig.add_traces(go.Scatter(  x=(start, end, end, start, start), y=(6,6,5,5,6),
                                        mode='markers+lines',
                                        line=dict(color=colors[features.index(row['feature'])]),
                                        name            =   'Annotations',
                                        hovertemplate   =   f'Coordinates: {end}-{start}' +
                                                            f'<br />Length: {end - start + 1}<br />' +
                                                            #f'<br />{"Strand"}: {row["strand"]}<br />' +
                                                            '<br />'.join([f'{x.title()}: {row[x]}' for x in ['strand', 'source', 'feature', 'score']]) +
                                                            f'<br />{"Attribute"}: {row["attribute"]}',
                                        fill            =   'toself',
                                        marker          =   dict(opacity=0),
                                        line_width      =   0,
                                        showlegend      =   False,
                            ), rows=1, cols=1
            )
    """    
    for i, row in alignments.df.iterrows():
        fig.add_traces(go.Scatter(  x=(start, end, end, start, start), y=(6,6,5,5,6),
                                    mode='markers+lines',
                                    line=dict(color=colors[features.index(row['feature'])]),
                                    name            =   'Reference',
                                    hovertemplate   =   f'Coordinates: {end}-{start}' +
                                                        f'<br />Length: {end - start + 1}<br />' +
                                                        #f'<br />{"Strand"}: {row["strand"]}<br />' +
                                                        '<br />'.join([f'{x.title()}: {row[x]}' for x in ['strand', 'source', 'feature', 'score']]) +
                                                        f'<br />{"Attribute"}: {row["attribute"]}',
                                    fill            =   'toself',
                                    marker          =   dict(opacity=0),
                                    line_width      =   0,
                                    showlegend      =   False,
                                    

                        ), rows=2, cols=1
        )
    """
    fig.update_traces(connectgaps=False)
    fig.update_layout(
        hovermode="x unified",
        #xaxis =   go.layout.XAxis(dict(title = 'Coordinates')),
        #yaxis =   go.layout.YAxis(dict(showticklabels=False)),
        #height =  600,
    )
    return fig


"""
Plot setting hooks
"""

def pan_drag(plot, element):
    """Set default drag mode to pan using a plot hook"""
    fig = plot.state
    fig['layout']['dragmode'] = 'pan'
    return fig

def plot_limits(plot, element, xmin=None, ymin=None, xmax=None, ymax=None):
    plot.handles['x_range'].min_interval = xmin
    plot.handles['y_range'].min_interval = ymin
    plot.handles['x_range'].max_interval = xmax
    plot.handles['y_range'].max_interval = ymax

def no_legend(plot, element):
    """Set default drag mode to pan using a plot hook"""
    fig = plot.state
    fig['show_legend'] = False
    return fig