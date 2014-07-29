import matplotlib.pyplot as plt
import collections
import numpy as np
import math

import aUtils
import aColors
#from bind_n_seq_pipeline import yield_kmers
import rbns_utils

def plot_stack(sp, bin_edges, enrichments, k, protein, scale='linear', thresh='2sig'):
    assert scale in ['linear', 'log']
    if thresh == '2sig':
        thresh = np.mean(enrichments) + 2 *np.std(enrichments)
    bottoms = []
    for bin_left, bin_right in aUtils.pairwise(bin_edges):
        color2height = collections.Counter()
        width = bin_right - bin_left
        for (kmeri, kmer), enrichment in zip(enumerate(rbns_utils.yield_kmers(k)), enrichments):
            if enrichment > bin_right or enrichment < bin_left:
                continue
            color = aColors.protein_colors(kmer, protein, enrichment >= thresh)
            color2height[color] += 1
        if scale == 'linear':
            bottom = 0.
            for color in aColors.ordered_colors:
                if not color2height[color]:
                    continue
                sp.bar(bin_left, color2height[color], width=width, facecolor=color, edgecolor='none', bottom=bottom)
                bottom += color2height[color]
        elif scale == 'log':
            bottom = 0.
            for color in aColors.ordered_colors:
                if not color2height[color]:
                    continue
                sp.bar(bin_left, math.log(color2height[color] + bottom + 1) - math.log(bottom + 1),
                  width=width, facecolor=color,
                  edgecolor='none', bottom=math.log(bottom + 1))
                bottom += color2height[color]
        bottoms.append(bottom)
    return max(bottoms)


if __name__ == '__main__':
    bin_edges = [0, 1, 2, 3, 4]
    enrichments = [0.5, 0.5, 2.5, 3.5]
    fig1 = plt.figure()
    sp = fig1.add_subplot(111)
    k = 1
    protein = 'test'
    plot_stack(sp, bin_edges, enrichments, k, protein, 'log')
    fig1.savefig('test_plot.png')
