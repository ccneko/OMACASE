#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220331


"""
Generate report for OM dataset
By Claire Chung
"""


import logging
from collections import OrderedDict
from os import makedirs
from pathlib import Path

# avoid multiple ray init with modin.pandas 
import pandas
from omacase.opmap import Opmaps, get_opmaps
from omacase.param_parser import *
from omacase import utils

logger = logging.getLogger(__name__)


def text_report(args):
    opmaps = get_opmaps(args=args, stat_only=True)
    if args.timestamp:
        timestamp = f'{utils.timestamp()}-'
    else:
        timestamp = ''
    if args.output_prefix is None:
        parent_dir  =   f'{Path(opmaps.filepath).parents[0]}/report'
        filename    =   opmaps.filename
    else:
        parent_dir  =   str(Path(args.output_prefix).parents[0])
        filename    =   str(Path(args.output_prefix).name)
    if not Path(parent_dir).is_dir():
        makedirs(parent_dir)
    output_filepath = f'{parent_dir}/{timestamp}{filename}-report.tsv'
    opmaps_summary_df(opmaps).to_csv(   output_filepath,
                                        sep='\t', index=None)
    logger.info(f'Generating dataset summary report: {output_filepath}')


def generate_tsv_report(opmaps: Opmaps, output_prefix, transformation=False):
    """Generate TSV report"""
    df = pandas.DataFrame([opmaps_summary_dict(opmaps)])
    if transformation:
        df = df.T
    output_filename = '.'.join([output_prefix, 'tsv'])
    df.to_csv(output_filename, sep='\t', index=None)
    return


def opmaps_summary_dict(opmaps: Opmaps=None):
    om_info_dict = OrderedDict()
    if opmaps is not None:
        om_info_dict['Report generation datetime']              =   utils.timestamp(format="%Y/%m/%d %H:%M:%S")
        om_info_dict['Filename']                                =   str(opmaps.filename)
        om_info_dict['File type']                               =   opmaps.file_format.upper()
        om_info_dict['File format version']                     =   opmaps.file_version
        """
        if opmaps.file_format == 'bnx':
            if opmaps.file_version == '1.3':
                om_info_dict['Platform']                        =  ' Saphyr'
        elif opmaps.file_format == 'cmap':
            if opmaps.file_version == '0.1':
                om_info_dict['Platform']                        =  ' Irys'
        """
        om_info_dict['Label motif']                             =   '; '.join(opmaps.motifs)
        logger.debug(f'Label motif: {om_info_dict["Label motif"].lower()}')
        if  'cttaag' in om_info_dict['Label motif'].lower() or\
            'dle1' in om_info_dict['Label motif'].lower():
            om_info_dict['Label motif']                         =   f'DLE1 ({om_info_dict["Label motif"]})'
            om_info_dict['Platform']                            =  'Bionano Saphyr'
        elif 'ctcgag' in om_info_dict['Label motif'].lower() or\
            'bsssi' in om_info_dict['Label motif'].lower():
            om_info_dict['Label motif']                         =   f'BSSSI ({om_info_dict["Label motif"]})'
        elif 'gctcttc' in om_info_dict['Label motif'].lower() or\
            'bspqi' in om_info_dict['Label motif'].lower():
            om_info_dict['Label motif']                         =   f'BSPQI ({om_info_dict["Label motif"]})'
        n50 = opmaps.get_n50()
        n90 = opmaps.get_n90()
        om_info_dict['Map number']                              =   f'{opmaps.mapnum:,}'
        om_info_dict['Total length']                            =   f'{int(opmaps.total_length):,}'
        om_info_dict['N50 length (bp)']                         =   f'{n50[0]:,}'
        om_info_dict['L50 value']                               =   f'{n50[1]:,}'
        om_info_dict['N90 length (bp)']                         =   f'{n90[0]:,}'
        om_info_dict['L90 value']                               =   f'{n90[1]:,}'
        om_info_dict['Max length (bp)']                         =   f'{max(opmaps.lengths):,}'
        om_info_dict['Min length (bp)']                         =   f'{min(opmaps.lengths):,}'
        # to add max len, min len, max SNR, min SNR
        om_info_dict['Label density / 100 kbp'] =   f'{opmaps.label_density:.3f}'
        logger.info(f'Label density / 100 kbp: {opmaps.label_density:.3f}')
        logger.info(f'Total length: {om_info_dict["Total length"]}')
        if opmaps.file_format == 'bnx':
            logger.info(f'Average QX11: {opmaps.qx11_avg:.3f}')
            om_info_dict['Mean of mol-average label SNR']       =   f'{opmaps.qx11_avg:.3f}'
            om_info_dict['SD of mol-average label SNR']         =   f'{opmaps.qx11_sd:.3f}'
            om_info_dict['Mean of mol-average label intensity'] =   f'{opmaps.qx12_avg:.3f}'
            om_info_dict['SD of mol-average label intensity']   =   f'{opmaps.qx12_sd:.3f}'
            om_info_dict['Mean of mol-average mol SNR']         =   f'{opmaps.mol_snr_avg:.3f}'
            om_info_dict['SD of mol-average mol SNR']           =   f'{opmaps.mol_snr_sd:.3f}'
            om_info_dict['Mean of mol-average mol intensity']   =   f'{opmaps.mol_intensity_avg:.3f}'
            om_info_dict['SD of mol-average mol intensity']     =   f'{opmaps.mol_intensity_sd:.3f}'
            """
            chip_info_list = opmaps.df['chip'].unique().to_list()
            om_info_dict['Chip(s)']             = ','.join([chip_info.split(',')[0].strip('SN_') \
                                                    for chip_info in chip_info_list])
            om_info_dict['Run ID(s)']           = ','.join([chip_info.split(',')[1].strip('Run_') \
                                                    for chip_info in chip_info_list])
            om_info_dict['No. of chip(s) used'] =   f"{len(om_info_dict['Chip(s)'])}"
            om_info_dict['No. of run(s)']       =   f"{len(om_info_dict['Run ID(s)'])}"
            #om_info_dict['Number of scans']         =   'unknown' # placeholder
            """
    return om_info_dict


def opmaps_summary_df(opmaps: Opmaps=None):
    om_info_dict    =   opmaps_summary_dict(opmaps)
    df              =   pandas.DataFrame(om_info_dict, index=[0]).T.reset_index()
    df.columns      =   ['Item', 'Dataset value']
    return df

