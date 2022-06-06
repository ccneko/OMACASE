#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220501


"""
Call tandem repeats from optical map data
By Claire Chung
"""


import logging

from intervaltree import IntervalTree
import numpy as np

from omacase.param_parser import *
from omacase.opmap import labels_to_segments, get_opmaps

logger = logging.getLogger(__name__)


def find_tandem_repeats_from_segments_ivtree(
                                segments, max_diff  =   MAX_DIFF,
                                max_diff_ratio      =   MAX_DIFF_RATIO,
                                max_avg_diff        =   MAX_AVG_DIFF,
                                max_avg_diff_ratio  =   MAX_AVG_DIFF_RATIO,
                                k                   =   K,
                                min_copy_num        =   MIN_COPY_NUM,
                                min_seg_len         =   MIN_SEG_LEN,
                                min_casette_length  =   MIN_CASETTE_LENGTH,
                                max_missed_label    =   MAX_MISSED_LABEL,
                                ivtree              =   None,
                                return_repeat_unit  =   False,
                            ):
    """Call tandem repeats from list of OM segment lengths"""
    copy_num = 1
    avg_diff = 0
    if return_repeat_unit:
        patterns = []

    for seg in range(k,len(segments)):
        offset = 0
        missed_label = 0
        if segments[seg] < min_seg_len:
            offset += segments[seg]
            missed_label   +=  1
            continue
        if sum(segments[seg-k-missed_label:seg]) < min_casette_length:
            continue
        diff = abs(segments[seg] - segments[seg - k - missed_label] - offset)
        curr_diff_ratio =   abs((float(sum(segments[seg-k-missed_label-1:seg-1]))\
                                -float(sum(segments[seg-k-missed_label:seg])))\
                                /float(sum(segments[seg-k-missed_label:seg])))
        if (diff <= max_diff) and (avg_diff <= max_avg_diff) \
            and curr_diff_ratio < max_diff_ratio \
            and missed_label <= max_missed_label:
            if seg!=0 and (seg - missed_label) % k == 0:
                offset = 0
                missed_label = 0
                copy_num += 1
            avg_diff = (avg_diff * (copy_num - 1) + diff)/float(copy_num)
        else:
            if copy_num >= min_copy_num:
                start_label = seg - k*(copy_num) - missed_label
                end_label = seg - missed_label
                
                repeat_unit_size = int(round(float(sum(segments[start_label:end_label]))/float(copy_num), 0))
                pattern = segments[start_label:end_label]
                if k==2:
                    if np.all([abs(segment - np.mean(pattern))< RESOLUTION/2 for segment in pattern]):
                        k = 1
                        repeat_unit_size /= 2
                        copy_num *= 2
                repeat_score = get_repeat_score(k, copy_num, segments=segments[start_label:end_label])
                if return_repeat_unit:
                    ivtree.addi(start_label, end_label, [k, copy_num, repeat_score, repeat_unit_size, pattern])
                else:
                    ivtree.addi(start_label, end_label, [k, copy_num, repeat_score, repeat_unit_size])
            copy_num = 1
            offset = 0
            missed_label = 0

    ivtree.merge_overlaps(strict=True, data_reducer=select_repeat_metrics)
    return ivtree


def find_tandem_repeats_from_opmaps_ivtree(
                                opmaps, max_diff    =   MAX_DIFF,
                                max_diff_ratio      =   MAX_DIFF_RATIO,
                                max_avg_diff        =   MAX_AVG_DIFF,
                                max_avg_diff_ratio  =   MAX_AVG_DIFF_RATIO,
                                k                   =   K,
                                min_copy_num        =   MIN_COPY_NUM,
                                min_seg_len         =   MIN_SEG_LEN,
                                min_casette_length  =   MIN_CASETTE_LENGTH,
                                max_missed_label    =   MAX_MISSED_LABEL,
                                threads             =   1,
                                return_repeat_unit  =   False,
                            ):
    """Call tandem repeats from Opmaps object"""
    repeat_dict = {}
    for row in opmaps.df.iterrows():
        omid = row[0]
        if omid not in repeat_dict:
            repeat_dict[omid] = IntervalTree()
        segments = labels_to_segments(row[1]['labels'], start=True)

        for j in range(2,k+1):
            results = find_tandem_repeats_from_segments_ivtree(
                                    segments            =   segments,
                                    max_diff            =   max_diff,
                                    max_diff_ratio      =   max_diff_ratio,
                                    max_avg_diff        =   max_avg_diff,
                                    max_avg_diff_ratio  =   max_avg_diff_ratio,
                                    k                   =   j,
                                    min_copy_num        =   min_copy_num,
                                    min_seg_len         =   min_seg_len,
                                    min_casette_length  =   min_casette_length,
                                    max_missed_label    =   max_missed_label,
                                    ivtree              =   repeat_dict[omid],
                                    return_repeat_unit  =   return_repeat_unit,
            )
        
            repeat_dict[omid] |= results
        repeat_dict[omid].merge_overlaps(strict=True, data_reducer=select_repeat_metrics)
        if len(repeat_dict[omid]) == 0:
            del repeat_dict[omid]
        else:
            logger.debug(f'Map ID: {omid}')
    logger.debug(f'repeat_dict: {repeat_dict}')
    return repeat_dict


def write_repeat_bed(args=None, opmaps=None, return_repeat_unit=True, output_prefix=None):
    if opmaps is None:
        opmaps = get_opmaps(args, labels_only=True)

    logger.info('Calling repeats')
    repeat_dict = find_tandem_repeats_from_opmaps_ivtree(opmaps, return_repeat_unit=return_repeat_unit)
    if output_prefix is not None:
        output_prefix += '-repeats'
    elif args is not None and hasattr(args, 'output_prefix') and args.output_prefix is not None:
        output_prefix = f'{args.output_prefix}-repeats'
    elif args is not None and hasattr(args, 'input_file') and args.input_file is not None:
        output_prefix = f'{args.input_file}-repeats'
    else:
        output_prefix = 'repeats'
        
    output_filepath = f'{output_prefix}.bed'
    with open(output_filepath, 'w') as f:
        logger.info(f'Writing repeat entries to {output_filepath}')
        for omid in repeat_dict.keys():
            for iv in repeat_dict[omid]:
                start = opmaps.df.loc[omid, "labels"][iv[0]-1]
                end = opmaps.df.loc[omid, "labels"][iv[1]-1]
                k, copy_num, score = iv[2][0], iv[2][1], iv[2][2]
                repeat_unit_size = iv[2][3]
                pattern = iv[2][4]
                f.write('\t'.join(map(str, [omid, start, end, iv[0], iv[1], k, copy_num, repeat_unit_size, f'{score:.2f}', pattern]))+'\n')
    return


def select_repeat_metrics(a,b):
    # take smaller k
    if a[0]<b[0]:
        return a
    elif b[0]<a[0]:
        return b
    # take higher repeat count
    elif a[0]==b[0]:
        if a[1]>b[1]:
            return a
        else:
            return b


def get_repeat_score(k=0, copy_num=0, segments=None):
    cov = np.std(segments)/np.mean(segments) * 100
    repeat_unit_size = sum(segments) / copy_num
    return k**2 + cov + np.log10(copy_num) + np.sqrt(repeat_unit_size/1000)

