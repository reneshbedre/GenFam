import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
from matplotlib import pyplot as plt
from anot.settings import BASE_DIR, MEDIA_ROOT
import os
import pandas as pd


def bar_enrich(_data_file, _user):
    data = pd.read_table(_data_file,
                         sep='\t',
                         header=0)

    # fam_name = data['Gene Family']
    fam_name = data['Short Name']
    pvalue = data['P-value']
    log10_pvalue = -(np.log10(data['P-value']))
    xaxis = np.arange(len(fam_name))
    width = 0.5

    fig, ax = plt.subplots()
    # ax.margins(x=0.5)
    # plt.grid(True)
    ax.bar(xaxis,
           log10_pvalue,
           width,
           alpha=0.5,
           color='b')
    # for ticks position
    ax.set_xticks(xaxis + .05)
    # ax.set_xlim(-1, 6)
    # plt.xlim(xmin=1)
    # ax.set_xticks(xaxis)
    ax.set_xticklabels(fam_name,
                       rotation="90",
                       fontsize="8"
                       )
    plt.xlim([min(xaxis) - 0.25, max(xaxis) + 0.5])
    ax.set_ylabel('-log10(P-Value)', fontweight="bold", fontsize="8")
    ax.set_xlabel("Gene Families", fontweight="bold", fontsize="8")
    # it was added to show the xticks labels as it were long
    plt.tight_layout()
    _fam_out_enrich_image_file = os.path.join(MEDIA_ROOT, "output_files", _user, "fam_enrich_out.png")
    # plt.figure(figsize=(3, 4))      # w, h
    # fig.set_size_inches(5, 8)
    # fig.set_size_inches(7.47, 7.73)
    pl.savefig(_fam_out_enrich_image_file, format='png', dpi=200, bbox_inches='tight')
