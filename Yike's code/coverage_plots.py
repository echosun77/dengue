# vim: fdm=indent
'''
author:     Fabio Zanini
date:       18/01/22
content:    Plot coverage from virus BAM files (3 kids, DWS).
'''
import os
import sys
import glob
import numpy as np
import pandas as pd
from Bio import SeqIO
import pysam
import matplotlib.pyplot as plt


data_fdn = '../../data/'


feature_named = {
    'capsid protein': 'C',
    'membrane glycoprotein': 'M',
    'envelope protein': 'E',
    'nonstructural protein NS1': 'NS1',
    'nonstructural protein NS2A': 'NS2A',
    'nonstructural protein NS2B': 'NS2B',
    'nonstructural protein NS3': 'NS3',
    'nonstructural protein NS4A': 'NS4A',
    'nonstructural protein NS4B': 'NS4B',
    #'2K peptide': '2K',
    'nonstructural protein NS5': 'NS5',
}


if __name__ == '__main__':

    patient_names = os.listdir(data_fdn+'virus_bam_files/')
    bam_fns = glob.glob(data_fdn+'virus_bam_files/*/DENV*.bam')

    ref_fn = f'{data_fdn}references/Zhiyuan_refs/DENV1_for_viral_reads_DWS.gb'
    ref = SeqIO.read(ref_fn, 'gb')

    # Take only first read for each UMI
    # + strand (reverse read), - strand (fwd read)
    coverage = np.zeros((2, len(ref)), np.int32)
    for bam_fn in bam_fns:
        umis = set()
        with pysam.AlignmentFile(bam_fn) as bamfile:
            for read in bamfile:
                if (not read.has_tag('CB')) or (not read.has_tag('UB')):
                    continue
                cell_barcode = read.get_tag('CB')
                umi = read.get_tag('UB')

                if (cell_barcode, umi) in umis:
                    continue
                umis.add((cell_barcode, umi))

                coverage[int(~read.is_reverse), read.get_reference_positions()] += 1

    fig, axs = plt.subplots(
            3, 1, figsize=(6.5, 4), sharex=True,
            gridspec_kw={'height_ratios': [5, 1.1, 5]},
            )

    # Genome annotation
    genome_positions = set()
    ax = axs[1]
    ax.set_axis_off()
    for fea in ref.features:
        if fea.type != 'mat_peptide':
            continue
        name = feature_named.get(fea.qualifiers['product'][0], None)
        if name is None:
            continue
        start, end = fea.location.nofuzzy_start, fea.location.nofuzzy_end
        genome_positions.add(start)
        genome_positions.add(end)
        ax.add_patch(plt.Rectangle(
            (start, 0), end - start, 1, fc='k', ec='k', lw=2, alpha=0.3,
        ))
        yt = 0.5 + int((end - start) < 1000) - 2 * int('B' in name)
        ax.text(0.5 * (start + end), yt, name, ha='center', va='center')
    ax.set_ylim(-0.2, 1.2)

    # Coverage
    x = np.arange(len(coverage[0]))
    colors = ['steelblue', 'orangered']
    for i, ax in enumerate(axs[::2]):
        ax.fill_between(
                x,
                np.zeros(len(coverage[0])),
                0.5 + coverage[i], color=colors[i], alpha=0.5)
    ymax = 1.1 * max(axs[0].get_ylim()[1], axs[1].get_ylim()[1])
    for i, ax in enumerate(axs[::2]):
        ax.set_ylim(10, ymax)
        ax.set_yscale('log')
        for pos in genome_positions:
            ax.axvline(pos, lw=1, color='k', alpha=0.2)
        #ax.grid(True)
    axs[0].set_ylabel('+ strand\ncoverage [reads]')
    axs[-1].set_ylabel('- strand\ncoverage [reads]')
    axs[-1].invert_yaxis()
    axs[-1].set_xlabel('Position in DENV genome [nucleotides]')
    axs[-1].set_xlim(0, len(ref))
    fig.tight_layout(h_pad=0)

    for ext in ['svg', 'pdf', 'png']:
        kwargs = {}
        if ext == 'png':
            kwargs['dpi'] = 300
        fig.savefig(
            f'../../figures/viral_reads/coverage_mirrorplot.{ext}',
            **kwargs,
        )

    #FIXME: buggy?
    if False:
        covdiff = coverage[0] - coverage[1]
        switches = [0] + list((np.diff(covdiff > 0)).nonzero()[0]) + [len(covdiff)]
        fig, ax = plt.subplots(figsize=(6, 3))
        for i in range(len(switches) - 1):
            start, end = switches[i], switches[i+1]
            xi = x[start: end]
            yi = covdiff[start: end]
            color = 'orangered' if yi[0] < 0 else 'steelblue'
            ax.fill_between(xi, np.zeros(len(xi)), yi, color=color, alpha=0.5)
        ax.set_yscale('symlog', linthresh=1000)
        ax.grid(True, axis='x')
        ax.set_xlabel('Position in DENV genome [nucleotides]')
        ax.set_xlim(0, 10500)
        fig.tight_layout()

    plt.ion(); plt.show()
