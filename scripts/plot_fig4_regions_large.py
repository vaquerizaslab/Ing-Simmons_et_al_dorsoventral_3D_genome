import fanc
import fanc.plotting as fancplot
import logging
import matplotlib
import os
# import sys
import seaborn

matplotlib.use('agg')
logging.basicConfig(level=logging.INFO)


def plot_region(name, region):
    output_file = os.path.join("figures", "figure_4_panels", name + ".pdf")
    logging.info("Working on %s", name)
    logging.info("Will write output to %s", output_file)

    gd7_nc14_hic = fanc.load(os.path.join("data", "hic", "merged", "gd7-nc14", "hic",
                               "gd7-nc14_5kb.hic"), mode="r")
    gd7_nc14_hic_plot = fancplot.HicPlot(gd7_nc14_hic, norm="log", vmin=1e-03, vmax=1e-01,
                           draw_minor_ticks=False, title="gd7", max_dist='250kb')
    gd7_nc14_diff = fanc.load(os.path.join("data", "hic", "merged", "gd7-nc14", "hic",
                               "diff_control-nc14_gd7-nc14_5kb.hic"), mode="r")
    gd7_nc14_diff_plot = fancplot.HicPlot(gd7_nc14_diff, norm="lin", colormap='bwr_r',
                              vmin=-0.01, vmax=0.01, draw_minor_ticks=False, max_dist='250kb')

    Tollrm910_nc14_hic = fanc.load(os.path.join("data", "hic", "merged", "Tollrm910-nc14", "hic",
                               "Tollrm910-nc14_5kb.hic"), mode="r")
    Tollrm910_nc14_hic_plot = fancplot.HicPlot(Tollrm910_nc14_hic, norm="log", vmin=1e-03, vmax=1e-01,
                           draw_minor_ticks=False, title="Tollrm910", max_dist='250kb')
    Tollrm910_nc14_diff = fanc.load(os.path.join("data", "hic", "merged", "Tollrm910-nc14", "hic",
                               "diff_control-nc14_Tollrm910-nc14_5kb.hic"), mode="r")
    Tollrm910_nc14_diff_plot = fancplot.HicPlot(Tollrm910_nc14_diff, norm="lin", colormap='bwr_r',
                              vmin=-0.01, vmax=0.01, draw_minor_ticks=False, max_dist='250kb')

    Toll10B_nc14_hic = fanc.load(os.path.join("data", "hic", "merged", "Toll10B-nc14", "hic",
                               "Toll10B-nc14_5kb.hic"), mode="r")
    Toll10B_nc14_hic_plot = fancplot.HicPlot(Toll10B_nc14_hic, norm="log", vmin=1e-03, vmax=1e-01,
                           draw_minor_ticks=False, title="Toll10B", max_dist='250kb')
    Toll10B_nc14_diff = fanc.load(os.path.join("data", "hic", "merged", "Toll10B-nc14", "hic",
                               "diff_control-nc14_Toll10B-nc14_5kb.hic"), mode="r")
    Toll10B_nc14_diff_plot = fancplot.HicPlot(Toll10B_nc14_diff, norm="lin", colormap='bwr_r',
                              vmin=-0.01, vmax=0.01, draw_minor_ticks=False, max_dist='250kb')

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

    h3k27me3_dict = {name: os.path.join("external_data", "koenecke_2016_2017", "chipseq_aligned",
                                       "H3K27me3_" + name + "_sorted_filtered_merged_canonical_chrs.bw")
                    for name in ["gd7", "tl10b"]}
    h3k27me3_dict["Tollrm910"] = os.path.join("external_data", "extra_chip-seq", "chipseq_aligned",
                                       "H3K27me3_Tollrm910_sorted_filtered_merged_canonical_chrs.bw")

    # ins_dict = {name: os.path.join("data", "boundaries", name + "_2kb_8.bw")
    #             for name in ["gd7-nc14", "Tollrm910-nc14", "Toll10B-nc14", "3-4h"]}

    rnaseq_ylim = fancplot.helpers.LimitGroup()
    rnaseq_ylim = [0, 10]
    h3k27ac_ylim = fancplot.helpers.LimitGroup()
    h3k27me3_ylim = fancplot.helpers.LimitGroup()

    rnaseq_plot_gd7 = fancplot.LinePlot(rnaseq_dict['gd7'], fill=False,
                         plot_kwargs={'color': "#648fff"},
                         draw_minor_ticks=False, aspect=0.05,
                         ylim=rnaseq_ylim, n_yticks=2)

    h3k27ac_plot_gd7 = fancplot.LinePlot(h3k27ac_dict['gd7'], fill=False,
                         plot_kwargs={'color': "#648fff"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27ac_ylim, n_yticks=2)
    h3k27me3_plot_gd7 = fancplot.LinePlot(h3k27me3_dict['gd7'], fill=False,
                         plot_kwargs={'color': "#648fff"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27me3_ylim, n_yticks=2)

    rnaseq_plot_Tollrm910 = fancplot.LinePlot(rnaseq_dict['tlrm910'], fill=False,
                         plot_kwargs={'color': "#dc267f"},
                         draw_minor_ticks=False, aspect=0.05,
                         ylim=rnaseq_ylim, n_yticks=2)

    h3k27ac_plot_Tollrm910 = fancplot.LinePlot(h3k27ac_dict['Tollrm910'], fill=False,
                         plot_kwargs={'color': "#dc267f"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27ac_ylim, n_yticks=2)
    h3k27me3_plot_Tollrm910 = fancplot.LinePlot(h3k27me3_dict['Tollrm910'], fill=False,
                         plot_kwargs={'color': "#dc267f"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27me3_ylim, n_yticks=2)

    rnaseq_plot_toll10b = fancplot.LinePlot(rnaseq_dict['tl10b'], fill=False,
                         plot_kwargs={'color': "#ffb000"},
                         draw_minor_ticks=False, aspect=0.05,
                         ylim=rnaseq_ylim, n_yticks=2)

    h3k27ac_plot_toll10b = fancplot.LinePlot(h3k27ac_dict['tl10b'], fill=False,
                          plot_kwargs={'color': "#ffb000"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27ac_ylim, n_yticks=2)
    h3k27me3_plot_toll10b = fancplot.LinePlot(h3k27me3_dict['tl10b'], fill=False,
                          plot_kwargs={'color': "#ffb000"},
                          draw_minor_ticks=False, aspect=0.05,
                          ylim=h3k27me3_ylim, n_yticks=2)

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

    plots = [gd7_nc14_hic_plot,
             gd7_nc14_diff_plot,
             rnaseq_plot_gd7,
             h3k27ac_plot_gd7, gd7_enh_plot,
             h3k27me3_plot_gd7,

             Tollrm910_nc14_hic_plot,
             Tollrm910_nc14_diff_plot,
             rnaseq_plot_Tollrm910,
             h3k27ac_plot_Tollrm910, Tollrm910_enh_plot,
             h3k27me3_plot_Tollrm910,

             Toll10B_nc14_hic_plot,
             Toll10B_nc14_diff_plot,
             rnaseq_plot_toll10b,
             h3k27ac_plot_toll10b, toll10b_enh_plot,
             h3k27me3_plot_toll10b,

             genes_plot
             ]

    with fancplot.GenomicFigure(plots, ticks_last=True) as gfig:
        fig, axes = gfig.plot(region)
        fig.savefig(output_file)


regions = [# example regions
           ("2L_example", "2L:5,500,000-6,000,000"),  # example region on 2L
           ("2R_example", "2R:20,500,000-21,000,000"), # drop in ssim in Toll10B
           ("3L_example", "3L:6,700,000-7,200,000"),
           ("3R_example", "3R:20,000,000-20,500,000"),
           ("3R_example2", "3R:22,250,000-23,250,000"),
           ("X_example", "X:9,500,000-10,000,000") # example region in X

           ]

for name, region in regions:
    plot_region(name, region)
