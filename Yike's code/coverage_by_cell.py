# vim: fdm=indent
'''
author:     Fabio Zanini
date:       18/01/22
content:    Plot coverage from virus BAM files (3 kids, DWS).
'''
import os
from pathlib import Path
import sys
import glob
from collections import Counter
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import seaborn as sns


data_fdn = '../../data/virus_bam_files/'


if __name__ == '__main__':

    patient_names = os.listdir(data_fdn)
    bam_fns = glob.glob(data_fdn+'*/*.bam')

    ratios_all = {}
    for bam_fn in bam_fns:
        patient = Path(bam_fn).parent.stem
        molecules = Counter()
        with pysam.AlignmentFile(bam_fn) as bamfile:
            for read in bamfile:
                if (not read.has_tag('CB')) or (not read.has_tag('UB')):
                    continue
                strand = '+' if read.is_reverse else '-'
                cell_barcode = read.get_tag('CB')
                umi = read.get_tag('UB')
                molecules[(cell_barcode, umi, strand)] += 1

        molecules = pd.Series(molecules)
        # Collapse UMI
        molecules[:] = 1
        # Conut UMI
        molecules = (molecules.reset_index()
                              .groupby(['level_0', 'level_2'])
                              .sum()
                              .unstack(1, fill_value=0)[0])

        # Find cells with largest #molecules
        cmost = molecules.sum(axis=1).nlargest(50)
        cmost = cmost.loc[cmost >= 8].index

        # Get ratio - / sum
        ratios = 1.0 * (molecules.loc[cmost].T / molecules.loc[cmost].sum(axis=1).values).T['-']

        ratios_all[patient] = ratios

    fig, ax = plt.subplots(figsize=(3, 3))
    colors = sns.color_palette('Set2', n_colors=3)
    for i, (patient, ratios) in enumerate(ratios_all.items()):
        x = list(np.sort(ratios.values))
        y = list(1.0 - np.linspace(0, 1, len(x)))
        x = [0] + x + [1]
        y = [1] + y + [0]
        ax.plot(x, y, label=patient, color=colors[i])
    ax.legend()
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.grid(True)
    ax.set_xlim(0, 1)
    ax.set_xlabel('- strand / vRNA within a cell')
    ax.set_ylabel('Cumulative over cells')
    ax.annotate(
        '',
        xy=(0.5, 0.2), xytext=(0.5, 0.5),
        arrowprops=dict(arrowstyle='-|>', lw=2),
        zorder=5,
        )
    ax.scatter([0.5], [0.5], s=80, color='k', zorder=5)

    fig.tight_layout()

    plt.ion(); plt.show()
