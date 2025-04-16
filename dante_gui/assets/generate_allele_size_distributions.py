#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from math import floor


def plot_allele_size_distribution(
    module_id, c, benign_r, premutation_r, pathogenic_r, motif_len, read_len=150
):
    counts = list(map(lambda x: int(x), c.split(sep=",")))
    benign_r = list(map(lambda x: int(x), benign_r.split("-")))
    premutation_r = list(map(lambda x: int(x), premutation_r.split("-")))
    pathogenic_r = list(map(lambda x: int(x), pathogenic_r.split("-")))

    benign_range = set()
    premutation_range = set()
    pathogenic_range = set()
    unknown_range = set(range(len(counts)))

    benign_range        |= set(range(benign_r[0], benign_r[1] + 1))
    premutation_range   |= set(range(premutation_r[0], premutation_r[1] + 1))
    pathogenic_range    |= set(range(pathogenic_r[0], pathogenic_r[1] + 1))
    unknown_range -= benign_range | premutation_range | pathogenic_range

    benign = [x if i in benign_range else 0 for i, x in enumerate(counts)]
    premutation = [x if i in premutation_range else 0 for i, x in enumerate(counts)]
    pathogenic = [x if i in pathogenic_range else 0 for i, x in enumerate(counts)]
    unknown = [x if i in unknown_range else 0 for i, x in enumerate(counts)]

    linewidth = 0.8
    cgreen = "#CEF1CA"
    cyellow = "#FFF1D3"
    cred = "#F9B4C8"
    cgray = "#D1D5E3"
    cdark = "#00000088"

    x = np.arange(len(counts))
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.bar(x, benign, color=cdark)
    ax.bar(x, premutation, color=cdark)
    ax.bar(x, pathogenic, color=cdark)
    ax.bar(x, unknown, color=cdark)

    start = -0.5
    ax.set_xlim(left=start)
    tmp1 = ax.add_patch(plt.Rectangle((start, 0), ax.get_xlim()[1], ax.get_ylim()[1], color=cgray, zorder=0))
    tmp2 = ax.add_patch(plt.Rectangle((start + benign_r[0], 0), benign_r[1] - benign_r[0] + 1, ax.get_ylim()[1], color=cgreen, zorder=0))
    tmp3 = ax.add_patch(plt.Rectangle((start + premutation_r[0], 0), premutation_r[1] - premutation_r[0] + 1, ax.get_ylim()[1], color=cyellow, zorder=0))
    tmp4 = ax.add_patch(plt.Rectangle((start + pathogenic_r[0], 0), pathogenic_r[1] - pathogenic_r[0] + 1 , ax.get_ylim()[1], color=cred, zorder=0))
    ax.legend([tmp2, tmp3, tmp4, tmp1], ["Benign", "Premutation", "Pathogenic", "Unknown"])

    vline_x = floor(read_len / motif_len) - start
    ax.axvline(vline_x, color="black", linewidth=linewidth)
    ax.text(vline_x, ax.get_ylim()[1] * 0.5, f"Read length ({read_len}bp)", rotation=90, horizontalalignment="right")
    ax.set_title(f"Allele size distribution for {module_id}")

    plt.show()
    plt.close()

# %%


module_id = "ALS-0"
c = \
    "2,0,0,0,1,2,4,6,7,9,"\
    "7,8,3,1,2,1,2,0,1,0,"\
    "1,2,4,2,3,3,1,2,0,0,"\
    "2,3,1,0,1,3,2,2,3,2,"\
    "1,1,0,0,0,2,3,2,0,2,"\
    "2,0,0,0,0,1,0,0,0,1,"\
    "1,0,0,0,0,2,0,1,0,2"

benign_r = "0-14"
premutation_r = "17-21"
pathogenic_r = "35-100"
motif_len = 5

plot_allele_size_distribution(module_id, c, benign_r, premutation_r, pathogenic_r, motif_len)
