import fanc
import fanc.plotting as fancplot
import logging
import matplotlib
import os
import sys
# import seaborn

from matplotlib.backends.backend_pdf import PdfPages

matplotlib.use('agg')
logging.basicConfig(level=logging.INFO)

output_prefix = sys.argv[1]

output_file = output_prefix + ".pdf"
logging.info("Will write output to %s", str(output_file))


def plot_regions(regions):

    h = fanc.load(os.path.join("data", "hic", "merged", "3-4h", "hic",
                               "3-4h_2kb.hic"), mode="r")
    h_plot = fancplot.HicPlot(h, vmin=1e-03, vmax=1e-01, norm="log",
                           draw_minor_ticks=False)

    genes = "external_data/flybase/dmel-all-r6.30.gtf.gz"
    genes_plot = fancplot.GenePlot(genes, squash=True, group_by="gene_symbol",
                                aspect=0.15, label_field="gene_symbol",
                                show_labels=False, draw_minor_ticks=False)

    rnaseq_dict = {name: os.path.join("external_data", "koenecke_2016_2017", "rnaseq_aligned",
                                      name + "_sorted_filtered_merged_canonical_chrs_rnaseq.bw")
                   for name in ["gd7", "tlrm910", "tl10b"]}
    h3k27ac_dict = {name: os.path.join("external_data", "koenecke_2016_2017", "chipseq_aligned",
                                       "H3K27ac_" + name + "_sorted_filtered_merged_canonical_chrs.bw")
                    for name in ["gd7", "tl10b"]}
    h3k27ac_dict["Tollrm910"] = os.path.join("external_data", "extra_chip-seq", "chipseq_aligned",
                   "H3K27ac_Tollrm910_sorted_filtered_merged_canonical_chrs.bw")

    rnaseq_ylim = fancplot.helpers.LimitGroup()
    h3k27ac_ylim = fancplot.helpers.LimitGroup()
    polii_ylim = fancplot.helpers.LimitGroup()

    # polii_chip_early = os.path.join("external_data", "blythe_2015", "aligned",
    #                           "PolII-pSer5_NC14-early_sorted_filtered_merged_canonical_chrs.bw")
    # polii_chip_mid = os.path.join("external_data", "blythe_2015", "aligned",
    #                           "PolII-pSer5_NC14-middle_sorted_filtered_merged_canonical_chrs.bw")
    polii_chip_late = os.path.join("external_data", "blythe_2015", "aligned",
                              "PolII-pSer5_NC14-late_sorted_filtered_merged_canonical_chrs.bw")

    # polii_early_plot = fancplot.LinePlot(polii_chip_early, fill=True, plot_kwargs={'color': "black"},
    #                                  draw_minor_ticks=False, aspect=0.05,
    #                                  ylim=polii_ylim, n_yticks=2)
    # polii_mid_plot = fancplot.LinePlot(polii_chip_mid, fill=True, plot_kwargs={'color': "black"},
    #                                  draw_minor_ticks=False, aspect=0.05,
    #                                  ylim=polii_ylim, n_yticks=2)
    polii_late_plot = fancplot.LinePlot(polii_chip_late, fill=True, plot_kwargs={'color': "black"},
                                     draw_minor_ticks=False, aspect=0.05,
                                     ylim=polii_ylim, n_yticks=2)

    rnaseq_plot_gd7 = fancplot.LinePlot(rnaseq_dict['gd7'], fill=True,
                         plot_kwargs={'color': "#648fff"},
                         draw_minor_ticks=False, aspect=0.05,
                         n_yticks=2)

    h3k27ac_plot_gd7 = fancplot.LinePlot(h3k27ac_dict['gd7'], fill=True,
                         plot_kwargs={'color': "#648fff"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27ac_ylim, n_yticks=2)

    rnaseq_plot_Tollrm910 = fancplot.LinePlot(rnaseq_dict['tlrm910'], fill=True,
                         plot_kwargs={'color': "#dc267f"},
                         draw_minor_ticks=False, aspect=0.05,
                         n_yticks=2)

    h3k27ac_plot_Tollrm910 = fancplot.LinePlot(h3k27ac_dict['Tollrm910'], fill=True,
                         plot_kwargs={'color': "#dc267f"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27ac_ylim, n_yticks=2)

    rnaseq_plot_toll10b = fancplot.LinePlot(rnaseq_dict['tl10b'], fill=True,
                         plot_kwargs={'color': "#ffb000"},
                         draw_minor_ticks=False, aspect=0.05,
                         n_yticks=2)

    h3k27ac_plot_toll10b = fancplot.LinePlot(h3k27ac_dict['tl10b'], fill=True,
                          plot_kwargs={'color': "#ffb000"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27ac_ylim, n_yticks=2)

    gd7_enh = "data/supplementary_tables/gd7_candidate_enhancers.bed"
    gd7_enh_plot = fancplot.GenomicFeaturePlot(gd7_enh,
                                            aspect=0.02, color="#648fff",
                                            draw_minor_ticks=False)

    Tollrm910_enh = "data/supplementary_tables/Tollrm910_candidate_enhancers.bed"
    Tollrm910_enh_plot = fancplot.GenomicFeaturePlot(Tollrm910_enh,
                                            aspect=0.02, color="#dc267f",
                                            draw_minor_ticks=False)

    toll10b_enh = "data/supplementary_tables/Toll10B_candidate_enhancers.bed"
    toll10b_enh_plot = fancplot.GenomicFeaturePlot(toll10b_enh,
                                                        aspect=0.02, color="#ffb000",
                                                        draw_minor_ticks=False)

    plots = [h_plot,
             # ins_plot,
             # boundaries_plot,
             genes_plot,
             # hk_plot,
             # polii_early_plot, polii_mid_plot,
             polii_late_plot,
             rnaseq_plot_gd7,
             rnaseq_plot_Tollrm910,
             rnaseq_plot_toll10b,
             h3k27ac_plot_gd7,
             h3k27ac_plot_Tollrm910,
             h3k27ac_plot_toll10b,
             gd7_enh_plot,
             Tollrm910_enh_plot,
             toll10b_enh_plot
             ]

    with PdfPages(output_file) as pdf:
        with fancplot.GenomicFigure(plots, ticks_last=True) as gfig:
            for name, region, rnaseq_ylim in regions:
                logging.info(region)
                fig, axes = gfig.plot(region)
                axes[3].set_ylim([0, rnaseq_ylim])
                axes[4].set_ylim([0, rnaseq_ylim])
                axes[5].set_ylim([0, rnaseq_ylim])
                pdf.savefig()


regions = [
           # expressed in Toll10B
           ("twi", "2R:22,900,000-23,100,000", 100),  # twi
           ("sna", "2L:15,400,000-15,600,000", 100),  # sna
           ("if", "X:16,700,000-16,850,000", 10),   # if
           ("NetA", "X:14,550,000-14,800,000", 20),   # NetA, NetB
           ("sog", "X:15,500,000-15,700,000", 20),   # sog
           # expressed in gd7
           ("Doc1", "3L:8,950,000-9,100,000", 10),    # Doc1, Doc2
           ("pnr", "3R:15,990,000-16,100,000", 10),  # pnr
           ("C15", "3R:21,450,000-21,550,000", 4),  # C15
           ]


plot_regions(regions)
