import fanc
import fanc.plotting as fancplot
import logging
import matplotlib
import os
# import sys
import seaborn

matplotlib.use('agg')
logging.basicConfig(level=logging.INFO)


def plot_region(name, region, promoter):
    output_file = os.path.join("figures", "figure_5_panels", name + ".pdf")
    logging.info("Working on %s", name)
    logging.info("Will write output to %s", output_file)

    gd7_stg10_hic = fanc.load(os.path.join("data", "hic", "merged", "gd7-stg10", "hic",
                               "gd7-stg10_2kb.hic"), mode="r")
    gd7_stg10_hic_plot = fancplot.HicPlot(gd7_stg10_hic, vmin=1e-03, vmax=1e-01, norm="log",
                           draw_minor_ticks=False, title="gd7-stg10")

    Tollrm910_stg10_hic = fanc.load(os.path.join("data", "hic", "merged", "Tollrm910-stg10", "hic",
                               "Tollrm910-stg10_2kb.hic"), mode="r")
    Tollrm910_stg10_hic_plot = fancplot.HicPlot(Tollrm910_stg10_hic, vmin=1e-03, vmax=1e-01, norm="log",
                           draw_minor_ticks=False, title="Tollrm910-stg10")

    Toll10B_stg10_hic = fanc.load(os.path.join("data", "hic", "merged", "Toll10B-stg10", "hic",
                               "Toll10B-stg10_2kb.hic"), mode="r")
    Toll10B_stg10_hic_plot = fancplot.HicPlot(Toll10B_stg10_hic, vmin=1e-03, vmax=1e-01, norm="log",
                           draw_minor_ticks=False, title="Toll10B-stg10")

    genes = "../external_data/flybase/dmel-all-r6.30.gtf.gz"
    genes_plot = fancplot.GenePlot(genes, squash=True, group_by="gene_symbol",
                                aspect=0.15, label_field="gene_symbol",
                                show_labels=False, draw_minor_ticks=False)

    rnaseq_dict = {name: os.path.join("../external_data", "koenecke_2016_2017", "rnaseq_aligned",
                                      name + "_sorted_filtered_merged_canonical_chrs_rnaseq.bw")
                   for name in ["gd7", "tlrm910", "tl10b"]}

    h3k27ac_dict = {name: os.path.join("../external_data", "koenecke_2016_2017", "chipseq_aligned",
                                       "H3K27ac_" + name + "_sorted_filtered_merged_canonical_chrs.bw")
                    for name in ["gd7", "tl10b"]}
    h3k27ac_dict["Tollrm910"] = os.path.join("../external_data", "extra_chip-seq", "chipseq_aligned",
                                       "H3K27ac_Tollrm910_sorted_filtered_merged_canonical_chrs.bw")

    ins_dict = {name: os.path.join("data", "boundaries", name + "_2kb_8.bw")
                for name in ["gd7-nc14", "Tollrm910-nc14", "Toll10B-nc14", "3-4h"]}

    rnaseq_ylim = fancplot.helpers.LimitGroup()
    h3k27ac_ylim = fancplot.helpers.LimitGroup()
    polii_ylim = fancplot.helpers.LimitGroup()

    # polii_chip_early = os.path.join("../external_data", "blythe_2015", "aligned",
    #                           "PolII-pSer5_NC14-early_sorted_filtered_merged_canonical_chrs.bw")
    # polii_chip_mid = os.path.join("../external_data", "blythe_2015", "aligned",
    #                           "PolII-pSer5_NC14-middle_sorted_filtered_merged_canonical_chrs.bw")
    polii_chip_late = os.path.join("../external_data", "blythe_2015", "aligned",
                                   "PolII-pSer5_NC14-late_sorted_filtered_merged_canonical_chrs.bw")

    # polii_early_plot = fancplot.LinePlot(polii_chip_early, fill=True, plot_kwargs={'color': "black"},
    #                                  draw_minor_ticks=False, aspect=0.05,
    #                                  ylim=polii_ylim, n_yticks=2)
    # polii_mid_plot = fancplot.LinePlot(polii_chip_mid, fill=True, plot_kwargs={'color': "black"},
    #                                  draw_minor_ticks=False, aspect=0.05,
    #                                  ylim=polii_ylim, n_yticks=2)
    polii_late_plot = fancplot.LinePlot(polii_chip_late, fill=False, plot_kwargs={'color': "black"},
                                     draw_minor_ticks=False, aspect=0.05,
                                     ylim=polii_ylim, n_yticks=2)

    rnaseq_plot_gd7 = fancplot.LinePlot(rnaseq_dict['gd7'], fill=False,
                         plot_kwargs={'color': "#648fff"},
                         draw_minor_ticks=False, aspect=0.05,
                         ylim=rnaseq_ylim, n_yticks=2)

    h3k27ac_plot_gd7 = fancplot.LinePlot(h3k27ac_dict['gd7'], fill=False,
                         plot_kwargs={'color': "#648fff"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27ac_ylim, n_yticks=2)

    rnaseq_plot_Tollrm910 = fancplot.LinePlot(rnaseq_dict['tlrm910'], fill=False,
                         plot_kwargs={'color': "#dc267f"},
                         draw_minor_ticks=False, aspect=0.05,
                         ylim=rnaseq_ylim, n_yticks=2)

    h3k27ac_plot_Tollrm910 = fancplot.LinePlot(h3k27ac_dict['Tollrm910'], fill=False,
                         plot_kwargs={'color': "#dc267f"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27ac_ylim, n_yticks=2)

    rnaseq_plot_toll10b = fancplot.LinePlot(rnaseq_dict['tl10b'], fill=False,
                         plot_kwargs={'color': "#ffb000"},
                         draw_minor_ticks=False, aspect=0.05,
                         ylim=rnaseq_ylim, n_yticks=2)

    h3k27ac_plot_toll10b = fancplot.LinePlot(h3k27ac_dict['tl10b'], fill=False,
                          plot_kwargs={'color': "#ffb000"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27ac_ylim, n_yticks=2)

    gd7_enh = "data/candidate_enhancers/h3k27ac_gd7_csaw_intersect.bed"
    gd7_enh_plot = fancplot.GenomicFeaturePlot(gd7_enh,
                                            aspect=0.02, color="#648fff",
                                            draw_minor_ticks=False)

    Tollrm910_enh = "data/candidate_enhancers/h3k27ac_Tollrm910_csaw_intersect.bed"
    Tollrm910_enh_plot = fancplot.GenomicFeaturePlot(Tollrm910_enh,
                                            aspect=0.02, color="#dc267f",
                                            draw_minor_ticks=False)

    toll10b_enh = "data/candidate_enhancers/h3k27ac_Toll10B_csaw_intersect.bed"
    toll10b_enh_plot = fancplot.GenomicFeaturePlot(toll10b_enh,
                                                        aspect=0.02, color="#ffb000",
                                                        draw_minor_ticks=False)

    ins_plot = fancplot.LinePlot(ins_dict, fill=False,
                              colors={"gd7-nc14": "#648fff80", "Tollrm910-nc14": "#dc267f",
                              "Toll10B-nc14": "#ffb00080", "3-4h": "#00000080"},
                              draw_minor_ticks=False, aspect=0.05,
                              n_yticks=2)
    boundaries = "data/boundaries/3-4h_final_boundaries.bed"
    boundaries_plot = fancplot.GenomicFeaturePlot(boundaries,
                                               aspect=0.02, color="black",
                                               draw_minor_ticks=False)
    housekeeping_promoters = "../external_data/tracks/housekeeping_promoters.bed"
    hk_plot = fancplot.GenomicFeaturePlot(housekeeping_promoters,
                                       aspect=0.02, color="red",
                                       draw_minor_ticks=False)

    plots = [gd7_stg10_hic_plot,
             Tollrm910_stg10_hic_plot,
             Toll10B_stg10_hic_plot,
             boundaries_plot,
             # ins_plot,
             genes_plot, 
             # hk_plot,
             # polii_early_plot, polii_mid_plot,
             polii_late_plot,
             rnaseq_plot_gd7, rnaseq_plot_Tollrm910, rnaseq_plot_toll10b,
             h3k27ac_plot_gd7, gd7_enh_plot,
             h3k27ac_plot_Tollrm910, Tollrm910_enh_plot,
             h3k27ac_plot_toll10b, toll10b_enh_plot
             ]

    with fancplot.GenomicFigure(plots, ticks_last=True) as gfig:
        fig, axes = gfig.plot(region)
        seaborn.despine(ax=axes[2], bottom=True, left=True)
        seaborn.despine(ax=axes[4], bottom=True)
        seaborn.despine(ax=axes[5], bottom=True)
        fig.savefig(output_file)


regions = [
           # expressed in Toll10B
           ("twi", "2R:22,900,000-23,100,000", "2R:22984374-22986373"),  # twi
           ("sna", "2L:15,400,000-15,600,000", "2L:15477261-15479260"),  # sna
           ("if", "X:16,700,000-16,850,000", "X:16782435-16784434"),   # if
           ("NetA", "X:14,550,000-14,800,000", "X:14652929-14654928"),   # NetA, NetB
           ("sog", "X:15,500,000-15,700,000", "X:15632435-15634434"),   # sog
           # expressed in gd7
           ("Doc1", "3L:8,950,000-9,100,000", "3L:9040375-9042374"),    # Doc1, Doc2
           ("pnr", "3R:15,990,000-16,100,000", "3R:16025079-16027078"),  # pnr
           ("C15", "3R:21,450,000-21,550,000", "3R:21498986-21500985"),  # C15
           ]

for name, region, promoter in regions:
    plot_region(name, region, promoter)
