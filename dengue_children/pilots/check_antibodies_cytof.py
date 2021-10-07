# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/09/21
content:    Check antibody specificity
'''
import os
import sys
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == '__main__':

    print('Load giant CSV file')
    fn_data = '../../data/csv/20k_dat.csv'
    df = pd.read_csv(fn_data)

    # Filter out Zika
    df = df.loc[df['status'] != 'Z'].copy()

    # Divide sick/healthy
    df['sick'] = df['status'].replace({
        'SD': 'sick', 'D': 'sick', 'DWS': 'sick', 'C': 'healthy'
    })

    gby = df.groupby('sick')
    ntots = gby.size()
    groups = ['healthy', 'sick']

    # I think the data is logged already
    aBs = ['DENV_NS1', 'DENV_NS3']

    # Single thresholds
    thresholds = np.linspace(0, 2.4, 50)
    res = pd.DataFrame(
        np.zeros((len(thresholds) - 1, 4)),
        columns=pd.MultiIndex.from_tuples([
            ('DENV_NS1', 'healthy'), ('DENV_NS1', 'sick'),
            ('DENV_NS3', 'healthy'), ('DENV_NS3', 'sick')]),
        )
    for gname in groups:
        dfi = gby.get_group(gname)
        # Subsample the sick
        if gname == 'sick':
            idx = np.random.choice(
                    np.arange(len(dfi)),
                    ntots['healthy'],
                    replace=False)
            dfi = dfi.iloc[idx]
        for aB in aBs:
            tmp = np.histogram(dfi[aB], bins=thresholds)[0]
            tmp = tmp[::-1].cumsum()[::-1]
            res[(aB, gname)] = tmp

    # Get specificity
    for aB in aBs:
        res[(aB, 'spec')] = 1.0 * res[(aB, 'sick')] / res[aB].sum(axis=1)

    fig, axs = plt.subplots(1, 2, figsize=(6, 3), sharex=True, sharey=True)
    xm = 0.5 * (thresholds[1:] + thresholds[:-1])
    for aB, ax in zip(aBs, axs):
        ax.set_title(aB)
        ax.plot(xm, 0.1 + res[(aB, 'sick')], lw=2, color='k', label='# sick cells')
        ax.set_ylim(bottom=0.1)
        ax.set_yscale('log')
        ax2 = ax.twinx()
        ax2.plot(xm, res[(aB, 'spec')], lw=2, color='tomato', label='Specificity')
        ax2.set_ylim(0, 1)
    fig.tight_layout()


    # Check the sum of all three antibodies
    df['DENV_NS1+NS3'] = df[['DENV_NS1', 'DENV_NS3']].sum(axis=1)
    df['DENV_ENV+NS1'] = df[['DENV_ENV', 'DENV_NS1']].sum(axis=1)
    df['DENV_ENV+NS3'] = df[['DENV_ENV', 'DENV_NS3']].sum(axis=1)
    df['DENV_sum3'] = df[['DENV_NS1', 'DENV_NS3', 'DENV_ENV']].sum(axis=1)
    res = pd.DataFrame(
        np.zeros((len(thresholds) - 1, 8)),
        columns=pd.MultiIndex.from_tuples([
            ('DENV_NS1+NS3', 'healthy'), ('DENV_NS1+NS3', 'sick'),
            ('DENV_ENV+NS1', 'healthy'), ('DENV_ENV+NS1', 'sick'),
            ('DENV_ENV+NS3', 'healthy'), ('DENV_ENV+NS3', 'sick'),
            ('DENV_sum3', 'healthy'), ('DENV_sum3', 'sick')]),
        )
    cols = ['DENV_NS1+NS3', 'DENV_ENV+NS1', 'DENV_ENV+NS3', 'DENV_sum3']
    for gname in groups:
        dfi = gby.get_group(gname)
        # Subsample the sick
        if gname == 'sick':
            idx = np.random.choice(
                    np.arange(len(dfi)),
                    ntots['healthy'],
                    replace=False)
            dfi = dfi.iloc[idx]
        for col in cols:
            tmp = np.histogram(dfi[col], bins=thresholds)[0]
            tmp = tmp[::-1].cumsum()[::-1]
            res[(col, gname)] = tmp

    # Get specificity
    for col in cols:
        res[(col, 'spec')] = 1.0 * res[(col, 'sick')] / res[col].sum(axis=1)

    fig, axs = plt.subplots(2, 2, figsize=(6, 6), sharex=True, sharey=True)
    axs = axs.ravel()
    xm = 0.5 * (thresholds[1:] + thresholds[:-1])
    for col, ax in zip(cols, axs):
        ax.set_title(col)
        ax.plot(xm, 0.1 + res[(col, 'sick')], lw=2, color='k', label='# sick cells')
        ax.set_ylim(bottom=0.1)
        ax.set_yscale('log')
        ax2 = ax.twinx()
        ax2.plot(xm, res[(col, 'spec')], lw=2, color='tomato', label='Specificity')
        ax2.set_ylim(0, 1)
    fig.tight_layout()

    fig, axs = plt.subplots(2, 2, figsize=(6, 6), sharex=True, sharey=True)
    axs = axs.ravel()
    for col, ax in zip(cols, axs):
        ax.set_title(col)
        y = 0.1 + res[(col, 'sick')]
        x = res[(col, 'spec')]
        ax.scatter(x, y, alpha=0.8, color='tomato')
        ax.set_ylim(bottom=0.1)
        ax.set_yscale('log')
        ax.set_xlim(0.5, 1)
        ax.grid(True)
    axs[3].set_xlabel('Specificity')
    axs[2].set_xlabel('Specificity')
    axs[2].set_ylabel('# sick cells')
    axs[0].set_ylabel('# sick cells')
    fig.tight_layout()

    # test a little bit of ML
    if False:
        # 1. balance
        idx = df['sick'] == 'healthy'
        tmp = np.random.choice(
                (df['sick'] == 'sick').values.nonzero()[0],
                idx.sum(),
                )
        idx.iloc[tmp] = True
        dfbal = df.loc[idx]

        # 2. get arrays
        X = dfbal[['DENV_NS1', 'DENV_NS3', 'DENV_ENV']].values
        y = (dfbal['sick'] == 'sick').astype(np.int16)

        # Split train/test
        from sklearn.model_selection import train_test_split, cross_val_score
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.33,
            )

        # 3. prepare model
        from sklearn.tree import DecisionTreeClassifier
        clf = DecisionTreeClassifier()

        # 4. run
        cvals = cross_val_score(
            clf, X, y, cv=5,
            fit_params=dict(sample_weight=10**dfbal['DENV_sum3'].values),
        )
        #clf.fit(X_train, y_train)

        # This doesn't really work because it's getting lost in those sick
        # cells without antibody staining...


    # double binning with some boolean gates
    from itertools import combinations
    res = {}
    pairs = list(combinations(['DENV_ENV', 'DENV_NS1', 'DENV_NS3'], 2))
    for gname in ['healthy', 'sick']:
        dfi = gby.get_group(gname)
        # Subsample the sick
        if gname == 'sick':
            idx = np.random.choice(
                    np.arange(len(dfi)),
                    ntots['healthy'],
                    replace=False)
            dfi = dfi.iloc[idx]
        for col1, col2 in pairs:
            resi = np.zeros((len(thresholds) - 1, len(thresholds) - 1), int)
            for i1, t1 in enumerate(thresholds[:-1]):
                for i2, t2 in enumerate(thresholds[:-1]):
                    # Antibody 1 is the rows (y axis!), antibody 2 is the columns (x axis!)
                    resi[i1, i2] = ((dfi[col1] >= t1) & (dfi[col2] >= t2)).sum()
            res[('# cells', gname, col1+'&'+col2)] = resi

    for col1, col2 in pairs:
        res[('# cells', 'total', col1+'&'+col2)] = \
            res[('# cells', 'healthy', col1+'&'+col2)] + \
            res[('# cells', 'sick', col1+'&'+col2)]
        res[('specificity', col1+'&'+col2)] = \
            res[('# cells', 'sick', col1+'&'+col2)] / \
            res[('# cells', 'total', col1+'&'+col2)]

    sensths = [300, 100, 30]
    fig, axs = plt.subplots(
            2 + len(sensths), 4, figsize=(10, 6 + 2 * len(sensths)),
            sharex=False, sharey=False,
            gridspec_kw=dict(width_ratios=[10] * 3 + [1]),
            )
    cmap = 'viridis'
    caxs = axs[:, -1]
    axs = axs[:, :-1]
    scm1 = plt.cm.ScalarMappable(plt.Normalize(vmin=0.5, vmax=1, clip=True), cmap=cmap)
    scm2 = plt.cm.ScalarMappable(plt.Normalize(vmin=-1, vmax=6, clip=True), cmap=cmap)
    for ip, (col1, col2) in enumerate(pairs):
        col = col1+'&'+col2
        axs[0, ip].set_title(col)
        sns.heatmap(
            res[('specificity', col)],
            cmap=cmap,
            vmin=0.5, vmax=1,
            cbar=False,
            ax=axs[0, ip]
        )
        sns.heatmap(
            np.log10(0.1 + res[('# cells', 'sick', col)]),
            vmin=-1, vmax=6,
            ax=axs[1, ip],
            cmap=cmap,
            cbar=False,
        )
        # Censored specificity
        for isens, sensth in enumerate(sensths):
            speccens = res[('specificity', col)].copy()
            speccens[res[('# cells', 'sick', col)] < sensth] = np.nan
            sns.heatmap(
                speccens,
                cmap=cmap,
                vmin=0.5, vmax=1,
                cbar=False,
                ax=axs[2+isens, ip]
            )
            # Max
            i1, i2 = np.unravel_index(np.nanargmax(speccens), speccens.shape)
            axs[2+isens, ip].scatter([i2], [i1], marker='*', color='red')
            axs[2+isens, ip].text(
                    i2+0.4, i1-0.2,
                    '{:.0%}\n{:}>={:.2f}\n{:}>={:.2f}\n#sick: {:}\n#healthy: {:}'.format(
                        speccens[i1, i2],
                        col1, thresholds[i1],
                        col2, thresholds[i2],
                        res[('# cells', 'sick', col)][i1, i2],
                        res[('# cells', 'healthy', col)][i1, i2],
                        ),
                    va='top',
                    color='red',
                    )

        for irow in range(2 + len(sensths)):
            axs[irow, ip].set_xticks([])
            axs[irow, ip].set_yticks([])
    cbar1 = fig.colorbar(scm1, cax=caxs[0])
    cbar2 = fig.colorbar(scm2, cax=caxs[1])
    cbar1.set_label('Specificity')
    cbar2.set_label('# cells from sick patients')
    cbar2.set_ticks([-1, 0, 1, 2, 3, 4, 5, 6])
    cbar2.set_ticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$', '$10^6$'])
    for isens, sensth in enumerate(sensths):
        cbar3 = fig.colorbar(scm1, cax=caxs[2 + isens])
        cbar3.set_label(f'Specificity [>= {sensth} sick cells]')
    fig.text(0.48, 0.02, 'Antibody 2', ha='center')
    fig.text(0.02, 0.48, 'Antibody 1', va='center', rotation=90)
    fig.tight_layout(rect=(0.04, 0.04, 1, 1))


