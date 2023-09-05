# vim: fdm=indent
'''
author:     Fabio Zanini
date:       18/01/22
content:    Plot coverage from human gene controls (3 kids, DWS)
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

sys.path.append('../../libraries')
import HTSeq


data_fdn = '../../data/virus_bam_files/'


def chromvector_to_array(cv):
     steps = list(cv.steps())
     if len(steps) < 3:
         return None, None
     start = steps[0][0].end
     array = np.zeros(steps[-1][0].start + 1 - start, int)
     for iv, count in steps[1:-1]:
         array[iv.start - start: iv.end - start] = count

     return array, start



if __name__ == '__main__':

     patient_names = os.listdir(data_fdn)
     bam_fns = glob.glob(data_fdn+'*/housekeep*_raw.bam')

     ## Find out what genes we have
     #for bam_fn in bam_fns:
     #    gene_counts = Counter()
     #    patient = Path(bam_fn).parent.stem
     #    with pysam.AlignmentFile(bam_fn) as bamfile:
     #        for ir, read in enumerate(bamfile):
     #            if (not read.has_tag('CB')) or (not read.has_tag('UB')):
     #                continue
     #            if read.has_tag('GN'):
     #                gene_counts[read.get_tag('GN')] += 1
     #    break

     #sys.exit()

     genes = [
         'ACTB',
         'TBP',
         'GAPDH',
         'ACTG1',
         'PSMD6',
         'PSMF1',
         ]
     # A few more that are less good (overlaps)


     ## Find out the coordinates
     #for bam_fn in bam_fns:
     #    patient = Path(bam_fn).parent.stem
     #    counts_by_gene = Counter()
     #    covs = {gene: HTSeq.GenomicArray('auto', typecode='i', stranded=True) for gene in genes}
     #    with HTSeq.BAM_Reader(bam_fn) as bamfile:
     #        for ir, read in enumerate(bamfile):
     #            #if (not read.has_tag('CB')) or (not read.has_tag('UB')):
     #            #    continue
     #            try:
     #                read.optional_field('CB')
     #                read.optional_field('UB')
     #                gene = read.optional_field('GN')
     #            except KeyError:
     #                continue

     #            if gene not in genes:
     #                continue

     #            if sum(counts_by_gene.values()) >= 1000 * len(genes):
     #                break

     #            if counts_by_gene[gene] >= 1000:
     #                continue

     #            counts_by_gene[gene] += 1
     #            covs[gene][read.iv] += 1
     #    break

     #fig, axs = plt.subplots(2, 3)
     #ygenes = {}
     #for gene, ax in zip(genes, axs.ravel()):
     #    for chrom, cd in covs[gene].chrom_vectors.items():
     #        yd = {}
     #        colors = {'+': 'steelblue', '-': 'darkred'}
     #        for strand, a in cd.items():
     #            y, start = chromvector_to_array(a)
     #            if y is None:
     #                continue
     #            yd[(chrom, strand)] = y, start
     #            x = start + np.arange(len(y))
     #            ax.plot(x, y, label=chrom+', '+strand, color=colors[strand])
     #        ax.set_ylim(bottom=0.9)
     #        ax.set_yscale('log')
     #        ax.legend()
     #    ax.set_title(gene)
     #    ygenes[gene] = yd
     #fig.tight_layout()

     #locations = {}
     #for gene in genes:
     #    yd = ygenes[gene]
     #    for (chrom, strand), (y, start) in yd.items():
     #        ind = (y > 0.1 * y.max()).nonzero()[0]
     #        loc = ind[0] + start, ind[-1] + 1 + start
     #        gene_strand = '-' if strand == '+' else '+'
     #        locations[gene] = (chrom, strand, loc[0], loc[1])
     #sys.exit()

     # Count reads within locations (no UMIs yet)
     gene_locations = {
         'ACTB': ('7', '+', 5527151, 5527254),
         'TBP': ('6', '-', 170554434, 170557021),
         'GAPDH': ('12', '-', 6534516, 6534859),
         'ACTG1': ('17', '+', 81509972, 81510107),
         'PSMD6': ('3', '+', 64010615, 64018670),
         'PSMF1': ('20', '-', 1118639, 1118783),
     }

     nreads_per_gene = 10000
     counts_by_gene_and_strand = Counter()
     for bam_fn in bam_fns:
         patient = Path(bam_fn).parent.stem
         print(patient)
         counts_by_gene = Counter()
         with pysam.AlignmentFile(bam_fn) as bamfile:
             for ir, read in enumerate(bamfile):
                 if (not read.has_tag('CB')) or (not read.has_tag('UB')):
                     continue

                 for gene, (chrom, gene_strand, start, end) in gene_locations.items():
                     if read.reference_name != chrom:
                         continue
                     tmp = read.get_reference_positions()
                     rstart, rend = tmp[0], tmp[-1] + 1
                     if rstart < start - 100:
                         continue
                     if rend > end + 100:
                         continue
                     break
                 else:
                     continue

                 if sum(counts_by_gene.values()) >= nreads_per_gene * len(genes):
                     break

                 if counts_by_gene[gene] >= nreads_per_gene:
                     continue

                 counts_by_gene[gene] += 1

                 rstrand = '-' if read.is_reverse else '+'
                 counts_by_gene_and_strand[(patient, gene, rstrand)] += 1

     counts = pd.Series(counts_by_gene_and_strand).unstack()
     fracs = (counts.T / counts.sum(axis=1)).T

     fracs = fracs.swaplevel(0, 1).sort_index()

     # Plot results
     genes_plot = ['ACTB', 'ACTG1', 'GAPDH', 'TBP', 'PSMF1']
     x = np.arange(len(genes_plot))
     fig, ax = plt.subplots(figsize=(3, 2))
     for i, gene in enumerate(genes_plot):
         y = fracs.loc[gene].values
         # Find minority
         y = y[:, y.mean(axis=0).argmin()]
         ax.scatter([i] * 3, y + 1e-5, color='k', alpha=0.5, s=60, zorder=5)
     ax.set_xticks(x)
     ax.set_xticklabels(genes_plot, rotation=90)
     ax.set_ylim(1e-4, 1e-1)
     ax.set_yscale('log')
     ax.grid(True)
     ax.set_ylabel('Fraction of reads\non opposite strand')
     fig.tight_layout()
     plt.ion(); plt.show()


     # Count reads within locations, with UMIs
     gene_locations = {
         'ACTB': ('7', '+', 5527151, 5527254),
         'TBP': ('6', '-', 170554434, 170557021),
         'GAPDH': ('12', '-', 6534516, 6534859),
         'ACTG1': ('17', '+', 81509972, 81510107),
         'PSMD6': ('3', '+', 64010615, 64018670),
         'PSMF1': ('20', '-', 1118639, 1118783),
     }

     nreads_per_gene = 10000
     counts_by_gene_and_strand = Counter()
     for bam_fn in bam_fns:
         patient = Path(bam_fn).parent.stem
         print(patient)
         with pysam.AlignmentFile(bam_fn) as bamfile:
             for ir, read in enumerate(bamfile):
                 if (not read.has_tag('CB')) or (not read.has_tag('UB')):
                     continue

                 for gene, (chrom, gene_strand, start, end) in gene_locations.items():
                     if read.reference_name != chrom:
                         continue
                     tmp = read.get_reference_positions()
                     rstart, rend = tmp[0], tmp[-1] + 1
                     if rstart < start - 100:
                         continue
                     if rend > end + 100:
                         continue
                     break
                 else:
                     continue

                 rstrand = '-' if read.is_reverse else '+'
                 cell_barcode = read.get_tag('CB')
                 umi = read.get_tag('UB')
                 counts_by_gene_and_strand[
                         (patient, gene, rstrand, cell_barcode, umi)
                     ] += 1

     counts = pd.Series(counts_by_gene_and_strand)
     # Collapse UMIs
     counts[:] = 1
     counts_umi = (counts.reset_index()
                         .groupby(['level_0', 'level_1', 'level_2'])
                         .sum()
                         .unstack()
                         .swaplevel(0, 1)
                         .sort_index()[0])

     fracs = (counts_umi.T / counts_umi.sum(axis=1)).T

     # Get viral reads
     counts_viral = Counter()
     viral_bam_fns = glob.glob(data_fdn+'*/DENV*.bam')
     for bam_fn in viral_bam_fns:
         patient = Path(bam_fn).parent.stem
         print(patient)
         with pysam.AlignmentFile(bam_fn) as bamfile:
             for ir, read in enumerate(bamfile):
                 if (not read.has_tag('CB')) or (not read.has_tag('UB')):
                     continue

                 rstrand = '-' if read.is_reverse else '+'
                 cell_barcode = read.get_tag('CB')
                 umi = read.get_tag('UB')
                 counts_viral[
                         (patient, rstrand, cell_barcode, umi)
                     ] += 1
     counts_viral = pd.Series(counts_viral)
     # Collapse UMI
     counts_viral[:] = 1
     # Count UMI
     counts_viral = (counts_viral.reset_index()
                           .groupby(['level_0', 'level_1'])
                           .sum()
                           .unstack()[0])
     fracs_viral = (counts_viral.T / counts_viral.sum(axis=1)).T

     # Plot results
     genes_plot = ['DENV', 'ACTB', 'ACTG1', 'GAPDH', 'TBP', 'PSMF1']
     x = np.arange(len(genes_plot))
     fig, ax = plt.subplots(figsize=(2.5, 2.3))
     for i, gene in enumerate(genes_plot):
         if gene != 'DENV':
             y = fracs.loc[gene].values
         else:
             y = fracs_viral.values
         # Find minority
         y = y[:, y.mean(axis=0).argmin()]

         color = 'tomato' if gene == 'DENV' else 'k'
         ax.scatter([i] * 3, y + 1e-5, color=color, alpha=0.5, s=60, zorder=5)

     ax.axhline(np.exp(np.log(fracs_viral.min(axis=1)).mean()), ls='--', lw=2, color='tomato')
     ax.axhline(np.exp(np.log(fracs.min(axis=1)).mean()), ls='--', lw=2, color='k')
     ax.set_xticks(x)
     ax.set_xticklabels(genes_plot, rotation=90)
     ax.set_ylim(1e-3, 1)
     ax.set_yscale('log')
     ax.grid(True)
     ax.set_ylabel('Fraction of reads\non opposite strand')
     fig.tight_layout()

     if False:
         for ext in ['svg', 'pdf', 'png']:
             kwargs = {}
             if ext == 'png':
                 kwargs['dpi'] = 300
             fig.savefig(
                 f'../../figures/viral_reads/fraction_molecules_opposite_strand.{ext}',
                 **kwargs,
             )

     plt.ion(); plt.show()