import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
from matplotlib import pyplot as plt
from anot.settings import BASE_DIR, MEDIA_ROOT
import os


def bar_enrich(_data_file, _user):
    # "/Users/renesh/Renesh_Docs/Research/django/anot/output_files/renesh/fam_enrich_out.txt"
    arr = np.genfromtxt(_data_file,
                        delimiter="\t",
                        dtype=None,
                        comments="#",
                        skip_header=1,
                        names=['fam', 'g1', 'g2', 'g3', 'g4', 'g5', 'pv', 'process', 'function', 'comp', 'fdr']
                        )
    fam_name = arr['fam']
    pvalue = arr['pv']
    log10_pvalue = -(np.log10(arr['pv']))
    xaxis = np.arange(len(fam_name))
    # print xaxis
    width = 0.25

    fig, ax = plt.subplots()
    plt.grid(True)
    ax.bar(xaxis,
           log10_pvalue,
           width,
           alpha=0.5,
           color='b')
    ax.set_xticks(xaxis + .1)
    # ax.set_xlim(-1, 6)
    # plt.xlim(xmin=1)
    # ax.set_xticks(xaxis)
    ax.set_xticklabels(fam_name,
                       rotation="90",
                       fontsize="8"
                       )
    plt.xlim([min(xaxis) - 0.25, max(xaxis) + 0.5])
    ax.set_ylabel('-log10(P-Value)', fontweight="bold", fontsize="16")
    ax.set_xlabel("Gene Families", fontweight="bold", fontsize="16")
    # it was added to show the xticks labels as it were long
    plt.tight_layout()
    _fam_out_enrich_image_file = os.path.join(MEDIA_ROOT, "output_files", _user, "fam_enrich_out.png")
    # plt.figure(figsize=(3, 4))      # w, h
    fig.set_size_inches(6, 5)
    pl.savefig(_fam_out_enrich_image_file, format='png', dpi=100)


# ar_enrich(data_file, user)