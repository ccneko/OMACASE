#!/usr/bin/env python3
# *_* coding: utf-8 *_*
# Last update: 20220516


"""
Simulate tandem repeats from input optical map (CMAP)
By Claire Chung
"""


import logging

from random import choice, choices, seed, shuffle, uniform
from scipy.stats import skewnorm
import numpy as np
import pandas as pd
from intervaltree import Interval, IntervalTree

from omacase.opmap import CMAP, labels_to_segments, get_opmaps
from omacase.param_parser import parse_opts_from_file, CN_LOC, CN_MAX, CN_MIN
from omacase.param_parser import CN_SCALE, CN_SKEW, K_LOC, K_MAX, K_MIN
from omacase.param_parser import K_SCALE, K_SKEW, LEN_ERR_MAX, RANDOM_SEED
from omacase.param_parser import REP_COUNT_LOC, REP_COUNT_MAX, REP_COUNT_MIN
from omacase.param_parser import REP_COUNT_SCALE, REP_COUNT_SKEW
from omacase.param_parser import REPEAT_UNIT_LENGTH_MIN, REPEAT_UNIT_LENGTH_MAX
from omacase.param_parser import SAMPLE_N, SEG_LENGTH_MAX, SEG_LENGTH_MIN
from omacase.param_parser import  TOTAL_REPEAT_LENGTH_MAX

logger = logging.getLogger(__name__)

repeats_cols = ['map_id', 'start', 'end', 'start_label', 'end_label', 'k', 'cn', 'score', 'repeat_unit_size', 'segments']

### Set the max number of repeats to simulate
sample_n = SAMPLE_N

### Set min and max of repeats to simulate in a map
rep_count_max = REP_COUNT_MAX
rep_count_min = REP_COUNT_MIN

### Set the per map repeat count skewed distribution
rep_count_skew = REP_COUNT_SKEW
rep_count_loc = REP_COUNT_LOC
rep_count_scale = REP_COUNT_SCALE

### Set min and max of number of segments in a repeat to simulate
k_min = K_MIN
k_max = K_MAX

### Set repeat unit segment number skewed distribution
k_skew = K_SKEW
k_loc = K_LOC
k_scale = K_SCALE

### Set min and max of copy number to simulate
cn_min = CN_MIN
cn_max = CN_MAX

### Set repeat unit segment number skewed distribution
cn_skew = CN_SKEW
cn_loc = CN_LOC
cn_scale = CN_SCALE

### Set segment length max error
len_err_max = LEN_ERR_MAX

### Set min and max repeat unit length
repeat_unit_length_min = REPEAT_UNIT_LENGTH_MIN
repeat_unit_length_max = REPEAT_UNIT_LENGTH_MAX

### Set min and max repeat unit fragment length (in bp)
seg_length_min = SEG_LENGTH_MIN
seg_length_max = SEG_LENGTH_MAX

### Set a maximum total repeat length (in bp) to add to each map
total_repeat_length_max = TOTAL_REPEAT_LENGTH_MAX

### Segments to exclude in each map (within repeats detected from the first round)
def repeat_sim(
                args=None, opmaps=None, included_source_maps=None, \
                output_prefix=None, random_seed=None, opt_file=None,
                excluded_regions_bed=None,
    ):
    """Optical map repeat simulation (CMAP)"""
    if opt_file is not None:
        args = parse_opts_from_file(opt_file)
    
    ## Set simulation parameters
    ### Set random seed for reproducibility
    if random_seed is not None:
        pass
    elif args is not None and args.seed is not None:
        random_seed = args.seed
    seed(random_seed)
    
    ## Set input and output files
    ### Load the input CMAP
    cmap = get_opmaps(args, labels_only=True)
    assert isinstance(cmap, CMAP)

    ### Set maps to include from choice of repeat simulation
    if included_source_maps is None:
        included_source_maps = cmap.df.index.to_list()

    ### Set segments to exclude from being drawn for repeat simulation
    to_exclude_trees = {}
    if args is not None and args.excluded_regions_bed is not None:
        excluded_regions_bed = args.excluded_regions_bed
    if excluded_regions_bed is not None:
        to_exclude_df = pd.read_csv(excluded_regions_bed, sep='\t', header=None, comment='#')[[0,1,2]]
        to_exclude_df.columns = ['seqname', 'start', 'end']
        for i, row in to_exclude_df.iterrows():
            if row['seqname'] not in to_exclude_trees:
                to_exclude_trees[row['seqname']] = IntervalTree()
            to_exclude_trees[row['seqname']].addi(row['start'], row['end'])

    ### Set output CMAP filename
    if output_prefix is not None:
        output_prefix += '-simulated_repeats'
    elif args is not None and hasattr(args, 'output_prefix') and args.output_prefix is not None:
        output_prefix = f'{args.output_prefix}-simulated_repeats'
    elif args is not None and hasattr(args, 'input_file') and args.input_file is not None:
        output_prefix = f'{args.input_file}-simulated_repeats'
    else:
        output_prefix = 'repeats'
    output_filepath = f'{output_prefix}.cmap'
    
    ### Simulation the max number of repeats per map
    rep_count_population = np.array([int(x) for x in skewnorm.rvs(rep_count_skew, loc=rep_count_loc, scale=rep_count_scale, size=2*len(included_source_maps), random_state=random_seed)])
    rep_count_population = rep_count_population[(rep_count_population >= rep_count_min) & (rep_count_population <= rep_count_max)]
    rep_count_population = choices(rep_count_population, k=len(included_source_maps))

    ### Simulate the repeat unit segment number (k) per repeat
    k_population = np.array([int(x) for x in skewnorm.rvs(k_skew, loc=k_loc, scale=k_scale, size=2*sample_n, random_state=random_seed)])
    k_population = k_population[(k_population >= k_min) & (k_population <= k_max)]
    k_population = choices(k_population, k=sample_n)

    ### Simulate the copy number (cn) per repeat
    cn_population = np.array([int(x) for x in skewnorm.rvs(cn_skew, loc=cn_loc, scale=cn_scale, size=2*sample_n, random_state=random_seed)])
    cn_population = cn_population[(cn_population >= cn_min) & (cn_population <= cn_max)]
    cn_population = choices(cn_population, k=sample_n)

    om_dict = {}
    has_repeat_map_dict = {}

    ## open an output file
    f = open(output_filepath, 'w')

    ## write header
    f.write('# hostname=omacase\n')
    f.write('# CMAP File Version:\t0.2\n')
    f.write('# Label Channels:\t1\n')
    f.write(f'# Number of Consensus Maps:\t{len(included_source_maps)}\n')
    f.write('#h\t' + '\t'.join(CMAP.default_columns) + '\n')
    f.write('#f\t' + '\t'.join(CMAP.col_dtype_desc) + '\n')

    dup_starts = {}
    dup_ends = {}

    shuffle(included_source_maps)

    ## Simulate the repeats
    total_repeats_simulated = 0
    for i, source_map in enumerate(included_source_maps):
        if total_repeats_simulated >= sample_n:
            rep_count_population[i:] = [0]*len(rep_count_population[i:])
        if rep_count_population[i] == 0:
            om_dict[source_map] = cmap.df.loc[source_map, 'labels']
        else:
            has_repeat_map_dict[source_map] = []
            labels = []
            source_label_pointer = 0
            per_map_repeats_simulated = 0
            #label_added = 0
            length_added = 0
            dup_starts[source_map] = []
            dup_ends[source_map] = []
            # 0-based source_label_pointer
            while per_map_repeats_simulated < rep_count_population[i] and source_label_pointer <  len(cmap.df.loc[source_map, 'labels']) - k_max - 1:
                ### Determine random start labels of simulated repeats
                k = k_population[total_repeats_simulated]
                cn = cn_population[total_repeats_simulated]
                source_label_pointers_candidates = range(int(round(source_label_pointer)), len(cmap.df.loc[source_map, 'labels']) - k_max)
                to_exclude_source_label_pointers = set()
                ### 0-based label_i
                for label_i, label in enumerate(cmap.df.loc[source_map, 'labels']):
                    if source_map in to_exclude_trees and to_exclude_trees[source_map].overlaps(label):
                        to_exclude_source_label_pointers.add(label_i)
                source_label_pointers_candidates = set(source_label_pointers_candidates) - to_exclude_source_label_pointers 
                logger.debug(source_label_pointers_candidates)
                ### 0-based chosen_start_source_label_pointer
                try:
                    chosen_start_source_label_pointer = choice(sorted(list(source_label_pointers_candidates)))
                except IndexError:
                    break
                logger.critical(f'source_map: {source_map}')
                logger.critical(f'chosen_start_source_label_pointer: {chosen_start_source_label_pointer}')
                logger.critical(cmap.df.loc[source_map, 'labels'][chosen_start_source_label_pointer])
                logger.debug(len(cmap.df.loc[source_map, 'labels'])-k_max)
                logger.debug(int(round(source_label_pointer)))
                logger.debug(cmap.df.loc[source_map, 'labels'][-1])
                ### 0-based dup_starts[source_map], dup_ends[source_map]
                dup_starts[source_map].append(chosen_start_source_label_pointer)
                dup_ends[source_map].append(dup_starts[source_map][-1] + k)
                
                ### Add labels until before the drawn repeat unit end, +1 to include the end label
                labels += [x + length_added for x in cmap.df.loc[source_map, 'labels'][source_label_pointer:dup_ends[source_map][-1]+1]]
                source_label_pointer = dup_ends[source_map][-1]
                ### Check if the segments meet simulation criteria 
                ### (roughly without considering the negligible max 5% error to be introduced for simplification)
                if np.any([x > seg_length_max for x in labels_to_segments(cmap.df.loc[source_map, 'labels'][dup_starts[source_map][-1]:dup_ends[source_map][-1]])]) or \
                    np.any([x < seg_length_min for x in labels_to_segments(cmap.df.loc[source_map, 'labels'][dup_starts[source_map][-1]:dup_ends[source_map][-1]])]):
                    logger.warning('Segment out of range. Skipped.')
                    logger.debug(labels_to_segments(cmap.df.loc[source_map, 'labels'][dup_starts[source_map][-1]-1:dup_ends[source_map][-1]]))
                    dup_starts[source_map].pop(-1)
                    dup_ends[source_map].pop(-1)
                    rep_count_population[i] -= 1
                    continue

                dup_start_coord = cmap.df.loc[source_map, 'labels'][dup_starts[source_map][-1]] + length_added
                logger.critical(f'labels: {labels[-1]}')
                source_label_pointer = dup_ends[source_map][-1]
                ### Simulate the repeat copies
                ### +1 to get k segments
                segments_to_repeat = labels_to_segments(cmap.df.loc[source_map, 'labels'][dup_starts[source_map][-1]:dup_ends[source_map][-1]+1])
                logger.critical(f'k: {k}')
                logger.critical(f'segments_to_repeat: {segments_to_repeat}')
                logger.critical(f'dup_starts[source_map][-1]: {dup_starts[source_map][-1]}')
                logger.critical(f'dup_ends[source_map][-1]+1: {dup_ends[source_map][-1]+1}')
                for repeat_copy_counter in range(cn-1):
                    for repeat_seg_pointer in range(k):
                        original_segment_length = segments_to_repeat[repeat_seg_pointer]
                        simulated_segment_length = round(original_segment_length * uniform(1-len_err_max, 1+len_err_max),1)
                        logger.critical(f'simulated_segment_length: {simulated_segment_length}')
                        length_added += simulated_segment_length
                        labels.append(labels[-1] + simulated_segment_length)
                source_label_pointer += k * cn
                dup_ends[source_map][-1] += k * (cn-1) ### added k before (line 200)
                per_map_repeats_simulated += 1
                total_repeats_simulated += 1
                size = round((round(labels[-1],1) - dup_start_coord)/cn, 0)
                ### print 1-based CMAP label index
                has_repeat_map_dict[source_map].append([dup_start_coord, round(labels[-1],1), dup_starts[source_map][-1]+1, dup_ends[source_map][-1]+1, k, cn, size])
                
            labels += [x + length_added for x in cmap.df.loc[source_map, 'labels'][source_label_pointer:]]
            if len(has_repeat_map_dict[source_map]) == 0:
                del has_repeat_map_dict[source_map]
            om_dict[source_map] = labels

    ## Print maps with simulated repeats to file
    for sim_map in sorted(om_dict.keys()):
        for i, label in enumerate(om_dict[sim_map]):
            label_num = len(om_dict[sim_map])
            map_length = f'{om_dict[sim_map][label_num-1]:.1f}'
            if sim_map in dup_starts and i in dup_starts[sim_map]:
                dup_start = 1
            else:
                dup_start = -1
            if sim_map in dup_ends and i in dup_ends[sim_map]:
                dup_end = 1
            else:
                dup_end = -1
            f.write('\t'.join([str (x) for x in [int(sim_map), map_length, label_num, i+1, 1, f'{label:.1f}'] \
                + [f'{-1:.2f}']*4 + [f'{dup_start:.2f}', f'{dup_end:.2f}'] + \
                [f'{0:.2f}']*3 + [f'{-1:.2f}', 0]])+'\n')
    f.close()
    sim_report_dicts = []
    for has_repeat_map in sorted(has_repeat_map_dict):
        for repeat in has_repeat_map_dict[has_repeat_map]:
            sim_report_dicts.append([has_repeat_map] + repeat)
    sim_df = pd.DataFrame(sim_report_dicts)
    sim_df.to_csv(f'{output_prefix}.bed', sep='\t', header=None, index=None)
    logger.info(f'Maps with simulated repeats written to {output_prefix}.cmap and {output_prefix}.bed')
    logger.info(f'Total number of simulated repeats: {total_repeats_simulated}')
    logger.info(f'Number of repeat-containing maps: {len(has_repeat_map_dict)}')
    logger.debug(f'Total number of input maps: {len(cmap.df.index)}')
    logger.info(f'Total number of output maps: {len(om_dict)}')
    return


def repeat_bed_to_ivtree(bed_filepath):
    repeat_df = pd.read_csv(bed_filepath, sep='\t', header=None)
    repeat_df.columns = repeats_cols
    repeat_trees = {}
    for map_id in repeat_df['map_id'].unique():
        repeat_trees[map_id] = IntervalTree()
        repeat_df.loc[repeat_df['map_id']==map_id].apply(lambda x: repeat_trees[map_id].addi(x['start'], x['end'], [x['k'], x['cn'], x['repeat_unit_size']]), axis=1)
    return repeat_trees