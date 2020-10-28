rule differential_h3k27ac_calculations:
    output:
        "scripts/differential_h3k27ac.html",
        expand("data/differential_chipseq/h3k27ac_{comparison}_{direction}.bed", 
            comparison=["gd7-tl10b", "gd7-Tollrm910", "tl10b-Tollrm910"],
            direction=["up", "down"])
    threads: 8
    shell:
        "R -e \"rmarkdown::render('/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/differential_H3K27ac.Rmd', "
        "output_file='/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/differential_H3K27ac.html', "
        "knit_root_dir='/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/')\""


rule make_enhancers:
    input:
        expand("data/differential_chipseq/h3k27ac_{comparison}_{direction}.bed", 
            comparison=["gd7-tl10b", "gd7-Tollrm910", "tl10b-Tollrm910"],
            direction=["up", "down"])
    output:
       "make_candidate_enhancers.html",
        expand("data/supplementary_tables/{sample}_candidate_enhancers.bed",
            sample=["gd7", "Tollrm910", "Toll10B"]),
        "figures/figure_1_panels/csaw_upset_plot.pdf",
        "figures/figure_1_panels/kvon_tile_activity.pdf"
    shell: 
        "R -e \"rmarkdown::render('/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/make_candidate_enhancers.Rmd', "
        "output_file='/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/make_candidate_enhancers.html', "
        "knit_root_dir='/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/')\""


rule make_differential_chipseq_hitile:
    input: 
        "data/supplementary_tables/candidate_enhancers/{track_id}.bed"
    output:
        "data/supplementary_tables/candidate_enhancers/{track_id}.hitile"
    shell: 
        "clodius aggregate bedfile "
        " --no-header "
        " --chromsizes-filename ../dm6_chrom_sizes_sanitized.txt "
        " -o {output} {input} "


# rule enhancer_domains_csaw:
#     input:
#         expand("data/candidate_enhancers/{track_id}.bed",
#             track_id = ["h3k27ac_gd7_csaw_intersect", "h3k27ac_Tollrm910_csaw_intersect", "h3k27ac_Toll10B_csaw_intersect"])
#     output:
#         expand("data/loops/enhancer_promoter_loops/{loop_set}.bedpe",
#          loop_set = ["gd7_csaw_ef_all", "gd7_csaw_ep_all", "Tollrm910_csaw_ef_all", 
#                 "Tollrm910_csaw_ep_all", "Toll10B_csaw_ef_all", "Toll10B_csaw_ep_all",
#                 "gd7_csaw_pf_all", "Tollrm910_csaw_pf_all", "Toll10B_csaw_pf_all"])
#     shell: 
#         "R -e \"rmarkdown::render('/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/csaw_enhancer_domains.Rmd', "
#         "output_file='/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/csaw_enhancer_domains.html', "
#         "knit_root_dir='/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/')\""


rule make_tf_peak_pairs:
    output:
        expand("data/TF_aggregates/{tf}_peak_pairs_{group}.bedpe",
            tf = ["Sna", "CBP", "Trl", "Pc", "Dl", "Mad", "Twi", "Zen"],
            group = ["all", "same_domain", "diff_domain", "100kb", "200kb", "500kb"])
    shell:
        "Rscript scripts/make_tf_peak_pairs.R"

rule make_tf_aggregates:
    input:
        mcool = "data/hic/cooler_files/{sample}_1kb.mcool",
        bedpe = "data/TF_aggregates/{tf}_peak_pairs_{group}.bedpe"
    output:
        "data/TF_aggregates/{sample}_1kb-1.0K_over_{tf}_peak_pairs_{group}_10-shifts.np.txt"
    threads: 8
    shell:
        "coolpup.py --n_proc {threads} --excl_chrs Y --outdir data/TF_aggregates/ " 
        "{input.mcool}::/resolutions/1000 {input.bedpe}"

rule plot_tf_aggregates:
    input:
        expand("data/TF_aggregates/{sample}_1kb-1.0K_over_{{tf}}_peak_pairs_{{group}}_10-shifts.np.txt",
            sample = ["control-nc14", "gd7-nc14", "Tollrm910-nc14", "Toll10B-nc14"])
    output:
        "figures/TF_aggregates/hi-c_all_1.0K_over_{tf}_peak_pairs_{group}_10-shifts.pdf"
    threads: 8
    shell:
        "plotpup.py --n_cols 1 "
        "--row_names control,gd7,Tollrm910,Toll10B "
        "--vmin 5e-1 --vmax 2e0 "
        "--enrichment 3 -o {output} {input} "


rule make_tf_aggregates_microc:
    input:
        mcool = "../micro-c/data/hic/cooler_files/{sample}_100bp.mcool",
        bedpe = "data/TF_aggregates/{tf}_peak_pairs_{group}.bedpe"
    output:
        "data/TF_aggregates/{sample}_100bp-1.6K_over_{tf}_peak_pairs_{group}_10-shifts.np.txt"
    threads: 8
    shell:
        "coolpup.py --n_proc {threads} --excl_chrs Y --outdir data/TF_aggregates/ " 
        "{input.mcool}::/resolutions/1600 {input.bedpe}"

rule plot_tf_aggregates_microc:
    input:
        expand("data/TF_aggregates/{sample}_100bp-1.6K_over_{{tf}}_peak_pairs_{{group}}_10-shifts.np.txt",
            sample = ["control_Rep1", "gd7_Rep1"])
    output:
        "figures/TF_aggregates/micro-c_all_1.6K_over_{tf}_peak_pairs_{group}_10-shifts.pdf"
    threads: 8
    shell:
        "plotpup.py --n_cols 1 "
        "--row_names control,gd7 "
        "--vmin 5e-1 --vmax 2e0 "
        "--enrichment 3 -o {output} {input} "


