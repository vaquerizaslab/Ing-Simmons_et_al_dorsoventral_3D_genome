import fanc
import logging
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
from fanc.plotting.statistics import saddle_plot

hic = snakemake.input['hic']
ab = snakemake.input['ab']
ev = snakemake.input['ev2']

plot = snakemake.output['plot']

n_percentiles = snakemake.params['n']
percentiles = [x for x in range(0, 101, 100 // int(n_percentiles))]

ab_matrix = fanc.load(ab, mode='r')
hic_matrix = fanc.load(hic, mode='r')
bed = fanc.load(ev)
ev = [r.score for r in bed.regions]

m, cutoffs = ab_matrix.enrichment_profile(hic_matrix,
                                          percentiles=percentiles,
                                          per_chromosome=True,
                                          only_gc=False,
                                          symmetric_at=None,
                                          exclude_chromosomes=(),
                                          eigenvector=ev,
                                          collapse_identical_breakpoints=True)


fig = plt.figure(figsize=(5, 5), dpi=300)
fig, axes = saddle_plot(m, cutoffs, colormap="RdBu_r", vmin=-1, vmax=1, only_gc=False,
                        fig=fig)
fig.savefig(plot)
plt.close(fig)
