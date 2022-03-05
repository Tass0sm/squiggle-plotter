#!/usr/bin/env python
# coding: utf-8

import os
import sys
import h5py
import pickle
from Bio import SeqIO
from poreplex import fast5_file
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import hashlib
import base64

# ?
from sklearn.linear_model import TheilSenRegressor
from scipy.interpolate import CubicSpline
import statsmodels
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.signal import medfilt
import random
import argparse

###############################################################################
#                                  Arguments                                  #
###############################################################################

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def file_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise NotADirectoryError(string)

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--reference',
                    type=file_path,
                    help='The reference fasta file.')

parser.add_argument('--sample-fast5-dir',
                    type=dir_path,
                    help='A directory containing fast5 files for the sample.')

parser.add_argument('--sample-scaling',
                    type=file_path,
                    help='A scaling parameter file for the sample.')

parser.add_argument('--control-fast5-dir',
                    type=dir_path,
                    help='A directory containing fast5 files for the control.')

parser.add_argument('--control-scaling',
                    type=file_path,
                    help='A scaling parameter file for the control.')

parser.add_argument('--target-position',
                    type=int,
                    help='The target position around which to center the output plot.')

parser.add_argument('--squiggle-width',
                    type=int,
                    help='The number of bases shown in the plot surrounding the target position.')

parser.add_argument('--output-file', '-o',
                    type=str,
                    help='The name of the output file (this can overwrite files).')

args = parser.parse_args()
print(args.output_file)

# Inputs:
# 1. The reference fasta
reference_fasta_path = args.reference

# 2. kmer-model: https://github.com/nanoporetech/kmer_models
kmer_model_path = "/users/PAS1405/tassosm/Desktop/kmer_models/r9.4_180mv_70bps_5mer_RNA/template_median69pA.model"

# 3. Scaling parameters. These were made by /users/PAS1405/tassosm/Desktop/squiggle-plot/fit-scaling-params.sh
scaling_ivt_path = args.control_scaling
scaling_vir_path = args.sample_scaling

# 4. The fast5 directories. I don't think these have to be from a guppy result directory.
ivt_fast5_dir = args.control_fast5_dir
vir_fast5_dir = args.sample_fast5_dir

###############################################################################
#                               setting up data                               #
###############################################################################

covseq = str(next(SeqIO.parse(open(reference_fasta_path), 'fasta')).seq)

kmer_model = pd.read_csv(kmer_model_path, sep='\t').set_index('kmer', drop=True)
kmer_mean = kmer_model['level_mean'].to_dict()

scaling_ivt = pd.read_csv(scaling_ivt_path, sep='\t', names=['read_id', 'scale', 'shift']).set_index('read_id')
scaling_vir = pd.read_csv(scaling_vir_path, sep='\t', names=['read_id', 'scale', 'shift']).set_index('read_id')
scaling = pd.concat([scaling_ivt, scaling_vir])

RESQUIGGLE_WIDTH = args.squiggle_width

def load_fast5_alignment(filename, start, end, corrected_group):
    read_id = filename.split('/')[-1].split('.')[0]

    with fast5_file.Fast5Reader(filename, read_id) as f5:
        h5 = f5.handle
        try:
            tnode = h5['Analyses/'+corrected_group+'/BaseCalled_template']
        except KeyError:
            return None

        if 'Alignment' not in tnode:
            return None;

        alnattrs = dict(tnode['Alignment'].attrs)

        if alnattrs['mapped_start'] > start or alnattrs['mapped_end'] < end:
            return None

        offset = tnode['Events'].attrs['read_start_rel_to_raw']
        rawsignal = f5.get_raw_data()[:-offset][::-1]
        events = pd.DataFrame(tnode['Events'][()])
        events['pos'] = events.index.to_series() + alnattrs['mapped_start']
        seq = b''.join(events['base']).decode()
        model_means = [kmer_mean[seq[left:left+5]] for left in range(len(seq)-4)]
        obs_means = events.iloc[2:-2].apply(
            lambda row: np.median(rawsignal[row['start']:row['start'] + row['length']]),
            axis=1)

        #regr = TheilSenRegressor(random_state=922)
        #regr.fit(np.array(obs_means)[:, np.newaxis], model_means)
        #scale, shift = regr.coef_[0], regr.intercept_
        scale, shift = scaling.loc[read_id][['scale', 'shift']]

        events = events[events['pos'].between(start, end - 1)]

        basepositions = []
        currentreadings = []

        for p, row in events.iterrows():
            sig = rawsignal[row['start']:row['start'] + row['length']] * scale + shift
            curpositions = np.linspace(row['pos'], row['pos'] + 1, len(sig) + 1)[:-1]
            basepositions.extend(curpositions)
            currentreadings.extend(sig)

        xplotpos = np.linspace(start, end, RESQUIGGLE_WIDTH * (end - start) + 1)
        csp = CubicSpline(basepositions, medfilt(currentreadings, 5))
        return xplotpos, csp(xplotpos)

# Make a dict from read ids (The hash in the name of a fast5 file) to the full fast5 file path.
ivt_fast5_map = {f.split('/')[-1].split('.')[0]: f for f in glob(os.path.join(ivt_fast5_dir, "*.fast5"))}
vir_fast5_map = {f.split('/')[-1].split('.')[0]: f for f in glob(os.path.join(vir_fast5_dir, "*.fast5"))}

###############################################################################
#                         Selecting Squiggle Location                         #
###############################################################################

# Edit for candidate positions.
# Candidate positions: 8079, 8975, 8989
TARGET_SITE = args.target_position
LEFT_END = TARGET_SITE - (args.squiggle_width // 2)
RIGHT_END = TARGET_SITE + (args.squiggle_width // 2)

target_read_ids = list(vir_fast5_map)
random.shuffle(target_read_ids)

control_read_ids = list(ivt_fast5_map)
random.shuffle(control_read_ids)

NUM_SIGNALS = 1000
CORRECTED_GROUP_IVT = 'RawGenomeCorrected_000'
CORRECTED_GROUP_VIRUS = 'RawGenomeCorrected_000'

virres_x = []; virres_y = []
ctlres_x = []; ctlres_y = []

def compute_xy(control_read_ids, target_read_ids):
    """Computes the x and y for the control signal and the target signal. Stores
    them in the global variables.
    """
    global ctlres_x, ctlres_y, virres_y, virres_x

    for read_id in target_read_ids:
        if read_id not in vir_fast5_map:
            continue
        filename = vir_fast5_map[read_id]
        res = load_fast5_alignment(filename, LEFT_END, RIGHT_END, CORRECTED_GROUP_VIRUS)
        if res is None:
            continue
        virres_x.append(res[0])
        virres_y.append(res[1])
        print(len(virres_x), end=' ')
        sys.stdout.flush()
        if len(virres_x) >= NUM_SIGNALS:
            break

    for read_id in control_read_ids:
        if read_id not in ivt_fast5_map:
            continue
        filename = ivt_fast5_map[read_id]
        res = load_fast5_alignment(filename, LEFT_END, RIGHT_END, CORRECTED_GROUP_IVT)
        if res is None:
            continue
        ctlres_x.append(res[0])
        ctlres_y.append(res[1])
        print(len(ctlres_x), end=' ')
        sys.stdout.flush()
        if len(ctlres_x) >= NUM_SIGNALS:
            break

def result_is_cached(pickle_name):
    return os.path.exists(pickle_name);

def compute_xy_memoized(control_read_ids, target_read_ids):
    """Tests if all inputs for this function have already been encountered. If so,
    fill the global variables with the cached data.
    """
    global ctlres_x, ctlres_y, virres_x, virres_y

    total_args = args.reference + args.control_fast5_dir + args.control_scaling + args.sample_fast5_dir + args.sample_scaling + str(args.target_position) + str(args.squiggle_width)
    total_hash_bytes = hashlib.sha224(bytes(total_args, "utf-8")).digest()
    total_hash_name = base64.b64encode(total_hash_bytes).decode('utf-8')
    signal_pickle_name = f"./saved-signals/{total_hash_name}-signals.pickle"

    if result_is_cached(signal_pickle_name):
        ctlres_x, ctlres_y, virres_x, virres_y = pickle.load(open(signal_pickle_name, 'rb'))
    else:
        compute_xy(control_read_ids, target_read_ids)
        pickle.dump([ctlres_x, ctlres_y, virres_x, virres_y], open(signal_pickle_name, 'wb'))

compute_xy_memoized(control_read_ids, target_read_ids)

###############################################################################
#                                    FIGURE                                   #
###############################################################################

XTICKS = np.arange(LEFT_END-1, RIGHT_END+1, 5)

fig, ax = plt.subplots(1, 1, figsize=(6.74, 3))

for i in range(35):
    if i == 0:
        ctlkwds = {'label': 'IVT'}; virkwds = {'label': 'Viral sgRNA S'}
    else:
        ctlkwds = virkwds = {}
    ax.plot(ctlres_x[i], ctlres_y[i], c='#000000', lw=.4, alpha=.4, zorder=3, **ctlkwds)
    ax.plot(virres_x[i], virres_y[i], c='#02a01e', lw=.4, alpha=.6, zorder=3, **virkwds)

for p in np.arange(LEFT_END, RIGHT_END):
    ax.axvline(p, c='black', alpha=.4, lw=.5, zorder=1)
    ax.annotate(covseq[p], (p + .5, 145), zorder=4, ha='center', va='top')

for sp in 'top right bottom'.split():
    ax.spines[sp].set_visible(False)

ax.set_xticks(XTICKS + .5)
ax.set_xticklabels(XTICKS + 1)
ax.xaxis.tick_top()
ax.set_ylabel('Normalized current (pA)')
ax.spines['left'].set_position(('outward', 5))
#plt.setp(ax.get_xticklines(), visible=False)

ax.set_ylim(75, 145)
#ax.legend(fontsize=10, loc='upper right')
ax.set_xlim(LEFT_END, RIGHT_END)
ax.grid(alpha=.2)
ax.spines['top'].set_position(('outward', 5))

plt.tight_layout()
plt.savefig(args.output_file)
