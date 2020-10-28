import os
import sys
import fanc
import fanc.plotting as fancplot
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pickle
import logging

logging.basicConfig(level=logging.INFO)

sample_name = sys.argv[1]

output_folder = "figures/figure_5_panels/"

logging.info("Working on sample %s", str(sample_name))

# set up data

hic_10kb_file = os.path.join("data", "micro-c", "merged", sample_name, "hic",
                             sample_name + "_10kb.hic")
hic_10kb = fanc.load(hic_10kb_file)

hic_1kb_file = os.path.join("data", "micro-c", "merged", sample_name, "hic",
                            sample_name + "_1kb.hic")
hic_1kb = fanc.load(hic_1kb_file)

# compartments data
hic_50kb_file = os.path.join("data", "micro-c", "merged", sample_name, "hic",
                             sample_name + "_50kb_masked.hic")

hic_50kb = fanc.load(hic_50kb_file)
ab_file = os.path.join("data", "micro-c", "merged", sample_name, "hic",
                       sample_name + "_50kb_masked.cor")
ab = fanc.ABCompartmentMatrix(ab_file, mode='r')
ab_enrichment_output_file = os.path.join(output_folder, sample_name + '_hic_50kb.ab_enrichment.npy')
ev_file = "data/compartments_by_gene_density/3-4h_50kb_masked_corrected_eigenvector.bed"
ev_regions = fanc.load(ev_file)
ev = [r.score for r in ev_regions.regions]
genome_file = "/home/research/vaquerizas/store/genomes/insects/Dmel/6.07/fasta/dmel-sanitized-chromosome-r6.07.fasta"

# tads and loops
tads_file = "data/boundaries/3-4h_paired_boundary_domains.bed"
tad_am_output_file = os.path.join(output_folder, sample_name + '_tad_aggregate.npy')

loops_file = "../external_data/cubenas-potts_2016/cubenas-potts_2016_loops.bedpe"
loop_am_output_file = os.path.join(output_folder, sample_name + '_loops_aggregate.npy')

large_plotting_region = fanc.GenomicRegion.from_string("2L:2,800,000-4,600,000")
small_plotting_region = fanc.GenomicRegion.from_string("2L:3,400,000-3,700,000")

# set up figure
matplotlib.rcParams.update({'font.size': 7})

fig = plt.figure(figsize=(9, 1.5), dpi=300)
gs = GridSpec(ncols=10, nrows=2,
              width_ratios=[10, 10, 10, 10, 2, 10, 1, 3, 1, 3],
              height_ratios=[8, 1],
              hspace=0.5, wspace=0.5)

# 1. plot large example region

ax_hic = plt.subplot(gs[0, 0])
p_hic = fancplot.HicPlot2D(hic_10kb, norm='log', vmin=1e-04, vmax=1e-01, ax=ax_hic, show_colorbar=False,
                           draw_tick_legend=False, draw_minor_ticks=False)
p_hic.plot(large_plotting_region)
ax_hic.set_xticks([large_plotting_region.start, large_plotting_region.end])
ax_hic.set_yticks([large_plotting_region.start, large_plotting_region.end])

# 2. plot compartment aggregate

if os.path.exists(ab_enrichment_output_file):
    with open(ab_enrichment_output_file, 'rb') as f:
        ab_enrichment_matrix, ev, cutoffs = pickle.load(f)
else:
    logging.info("Calculating compartment aggregate")
    ab_enrichment_matrix, cutoffs = ab.enrichment_profile(hic_50kb,
                                              percentiles=list(range(0, 100, 2)),
                                              per_chromosome=True,
                                              only_gc=False,
                                              symmetric_at=None,
                                              exclude_chromosomes=["3L", "3R"],
                                              eigenvector=ev,
                                              collapse_identical_breakpoints=True)

    with open(ab_enrichment_output_file, 'wb') as o:
        pickle.dump([ab_enrichment_matrix, ev, cutoffs], o)


ax_ab_enrichment = plt.subplot(gs[0, 1])
im_ab = ax_ab_enrichment.imshow(ab_enrichment_matrix, cmap='RdBu_r', vmin=-1, vmax=1,
                                interpolation='nearest', aspect='auto')

ax_ab_enrichment.set_xticks([0, ab_enrichment_matrix.shape[1] - 1])
ax_ab_enrichment.set_xticklabels(['A', 'B'])
# xlabels = ax_ab_enrichment.get_xticklabels()
# xlabels[0].set_horizontalalignment('left')
# xlabels[1].set_horizontalalignment('right')

ax_ab_enrichment.set_yticks([0, ab_enrichment_matrix.shape[1] - 1])
ax_ab_enrichment.set_yticklabels(['A', 'B'])
# ylabels = ax_ab_enrichment.get_yticklabels()
# ylabels[0].set_verticalalignment('top')
# ylabels[1].set_verticalalignment('bottom')


# 3. plot domain aggregate

tads = fanc.load(tads_file)
if not os.path.exists(tad_am_output_file):
    logging.info("Calculating TAD aggregate")
    tad_am = fanc.AggregateMatrix.from_regions(hic_1kb, tads.regions, log=True, oe=True,
                                               file_name=tad_am_output_file)
else:
    tad_am = fanc.AggregateMatrix(tad_am_output_file)

ax_aggregate_tads = plt.subplot(gs[0, 2])
tam = tad_am.matrix()
im_ab = ax_aggregate_tads.imshow(tam, cmap='RdBu_r', vmin=-1, vmax=1,
                                 interpolation='nearest', aspect='auto')

ylim = ax_aggregate_tads.get_ylim()
ax_aggregate_tads.set_xticks([tam.shape[0]/3, tam.shape[0]*2/3])
ax_aggregate_tads.set_yticks([tam.shape[0]/3, tam.shape[0]*2/3])
ax_aggregate_tads.set_xticklabels([])
ax_aggregate_tads.set_yticklabels([])
ax_aggregate_tads.set_ylim(ylim)

# 4. plot loop aggregate

loops = fanc.load(loops_file)
if not os.path.exists(loop_am_output_file):
    logging.info("Calculating loop aggregate")
    loop_am = fanc.AggregateMatrix.from_center_pairs(hic_1kb, loops, log=True, oe=True,
                                                     file_name=loop_am_output_file)
else:
    loop_am = fanc.AggregateMatrix(loop_am_output_file)

ax_aggregate_loops = plt.subplot(gs[0, 3])

lam = loop_am.matrix()
im_ab = ax_aggregate_loops.imshow(lam, cmap='RdBu_r', vmin=-1, vmax=1,
                                  interpolation='nearest', aspect='auto')

ylim = ax_aggregate_loops.get_ylim()
ax_aggregate_loops.set_xticks([lam.shape[0]/2])
ax_aggregate_loops.set_yticks([lam.shape[0]/2])
ax_aggregate_loops.set_xticklabels(['loop anchor'])
ax_aggregate_loops.set_yticklabels([''])
ax_aggregate_loops.set_ylim(ylim)

cax_aggregate_loops = plt.subplot(gs[0, 8])
cb_loops = plt.colorbar(im_ab, cax=cax_aggregate_loops)
cb_loops.set_label('Average log2-O/E\ncontacts')
cb_loops.set_ticks([-1, 0, 1])

# 5. plot smaller example region

ax_hic_small = plt.subplot(gs[0, 5])
cax_hic_small = plt.subplot(gs[0, 6])
p_hic_small = fancplot.HicPlot2D(hic_1kb, norm='log', vmin=1e-04, vmax=1e-01,
                                 ax=ax_hic_small, cax=cax_hic_small,
                                 draw_tick_legend=False, draw_minor_ticks=False)
p_hic_small.plot(small_plotting_region)
ax_hic_small.set_xticks([small_plotting_region.start, small_plotting_region.end])
ax_hic_small.set_yticks([small_plotting_region.start, small_plotting_region.end])

p_hic_small.add_colorbar(ax=cax_hic_small, aspect=40)
p_hic_small.colorbar.ax.minorticks_off()
p_hic_small.colorbar.set_ticks([1e-04, 1e-03, 1e-02, 1e-01])
p_hic_small.colorbar.set_label("Normalised contact\nprobability")


# save figure!

fig.savefig(os.path.join(output_folder, sample_name + '_overview_figure.pdf'))
