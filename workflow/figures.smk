###################
# FIGURES
###################

rule get_de_gene_regions:
	output:
		"external_data/koenecke_2016_2017/rnaseq_results/gd7_vs_tl10b_locations_of_de_genes_with_enhancers.bed",
		"external_data/koenecke_2016_2017/rnaseq_results/gd7_vs_tlrm910_locations_of_de_genes_with_enhancers.bed",
		"external_data/koenecke_2016_2017/rnaseq_results/tlrm910_vs_tl10b_locations_of_de_genes_with_enhancers.bed"
	script:
		"scripts/get_de_gene_regions.R"

rule plot_de_gene_regions_nc14:
	input:
		"external_data/koenecke_2016_2017/rnaseq_results/{rnaseq_comparison}_locations_of_de_genes_with_enhancers.bed",
	output:
		"figures/de_genes_regions/{rnaseq_comparison}_de_gene_regions_nc14.pdf"
	shell:
		"python scripts/plot_de_gene_regions.py {input} {output}"

rule plot_de_gene_regions_stg10:
	input:
		"external_data/koenecke_2016_2017/rnaseq_results/{rnaseq_comparison}_locations_of_de_genes_with_enhancers.bed",
	output:
		"figures/de_genes_regions/{rnaseq_comparison}_de_gene_regions_stg10.pdf"
	shell:
		"python scripts/plot_de_gene_regions_stg10.py {input} {output}"

## Figure 1

rule plot_fig1_regions:
	input:
		"data/hic/merged/3-4h/hic/3-4h_2kb.hic"
	output: 
		"figures/figure_1_panels/figure1_regions.pdf"
	params:
		output_prefix = "figures/figure_1_panels/figure1_regions"
	shell:
		"python scripts/plot_fig1_regions.py {params.output_prefix}"

rule enhancers_figures:
    input:
        expand("data/supplementary_tables/{sample}_candidate_enhancers.bed",
            sample=["gd7", "Tollrm910", "Toll10B"])
    output:
        "figures/figure_1_panels/csaw_gene_density_of_domains_with_enhancers.pdf",
        "figures/figure_1_panels/csaw_upset_plot.pdf",
        "figures/figure_1_panels/kvon_tile_activity.pdf",
        expand("data/loops/enhancer_promoter_loops/{loop_set}.bedpe",
         loop_set = ["gd7_csaw_ef_all", "gd7_csaw_ep_all", "Tollrm910_csaw_ef_all", 
                "Tollrm910_csaw_ep_all", "Toll10B_csaw_ef_all", "Toll10B_csaw_ep_all",
                "gd7_csaw_pf_all", "Tollrm910_csaw_pf_all", "Toll10B_csaw_pf_all"]),
        expand("data/loops/enhancer_enhancer_loops/{set}_enhancer_pairs_same_domain.bedpe",
            set=["gd7", "Tollrm910", "Toll10B"])
    shell:  
        "R -e \"rmarkdown::render('/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/plot_enhancers_figures.Rmd', "
        "output_file='/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/plot_enhancers_figures.html', "
        "knit_root_dir='/home/research/vaquerizas/liz/dorsal_ventral/for_paper/scripts/')\""

rule make_enhancer_heatmap_matrix_5kb:
    input: 
        expand("data/supplementary_tables/{sample}_candidate_enhancers.bed",
            sample=["gd7", "Tollrm910", "Toll10B"])
    output: "data/deeptools_matrices/candidate_enhancers_histone_mods_chipseq_5kb.gz"
    threads: 16
    shell:
        """
        computeMatrix reference-point --referencePoint center \
            --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
            --scoreFileName ../external_data/koenecke_2016_2017/chipseq_aligned/H3K27ac_gd7_sorted_filtered_merged_canonical_chrs.bw \
            external_data/extra_chip-seq/chipseq_aligned/H3K27ac_Tollrm910_sorted_filtered_merged_canonical_chrs.bw \
            external_data/koenecke_2016_2017/chipseq_aligned/H3K27ac_tl10b_sorted_filtered_merged_canonical_chrs.bw \
            external_data/koenecke_2016_2017/chipseq_aligned/H3K27me3_gd7_sorted_filtered_merged_canonical_chrs.bw \
            external_data/extra_chip-seq/chipseq_aligned/H3K27me3_Tollrm910_sorted_filtered_merged_canonical_chrs.bw \
            external_data/koenecke_2016_2017/chipseq_aligned/H3K27me3_tl10b_sorted_filtered_merged_canonical_chrs.bw \
            --samplesLabel H3K27ac_gd7 H3K27ac_Tollrm910 H3K27ac_Toll10B H3K27me3_gd7 H3K27me3_Tollrm910 H3K27me3_Toll10B \
            --regionsFileName {input} \
            --sortRegions keep \
            --outFileName {output} \
            -p {threads}
        """

rule make_enhancer_heatmap_5kb:
    input: "data/deeptools_matrices/candidate_enhancers_histone_mods_chipseq_5kb.gz"
    output: "figures/figure_1_panels/candidate_enhancers_histone_mods_chipseq_5kb.pdf"
    shell: 
        """
        plotHeatmap --matrixFile {input} \
        --outFileName {output} \
        --regionsLabel gd7 Tollrm910 Toll10B \
        --missingDataColor 1 \
        --colorMap Blues \
        --zMin 0 \
        --refPointLabel Enh \
        --whatToShow 'heatmap and colorbar' \
        --zMax 9 \
        --xAxisLabel ''\
        --heatmapHeight 14 \
        --heatmapWidth 2
        """

# Figure 4

rule plot_fig4_regions:
	input:
		"data/hic/merged/3-4h/hic/3-4h_2kb.hic"
	output: 
		expand("figures/figure_4_panels/{name}.pdf", 
			name=["twi", "sna", "if", "NetA", "sog", "Doc1", "pnr", "C15"])
	shell:
		"python scripts/plot_fig4_regions.py"

rule plot_fig4_regions_large:
	input:
		"data/hic/merged/3-4h/hic/3-4h_5kb.hic",
		expand("data/hic/merged/{sample}/hic/diff_control-nc14_{sample}_5kb.hic",
			sample=["gd7-nc14", "Tollrm910-nc14", "Toll10B-nc14"])
	output: 
		expand("figures/figure_4_panels/{name}.pdf", 
			name=["2L_example", "2R_example", "3L_example", "3R_example", "3R_example2", "X_example"])
	shell:
		"python scripts/plot_fig4_regions_large.py"

# Figure 3

rule plot_overview_figure:
	input: 
		ancient("data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_10kb.hic")
	output:
		"figures/figure_3_panels/{merged_sample_name}_overview_figure.pdf"
	shell:
		"python scripts/plot_overview_figure.py {wildcards.merged_sample_name}"

# Figure S4

rule plot_whole_genome:
	input:
		"data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}_wholegenomebalanced.hic"
	output:
		"figures/whole_genome/{merged_sample_name}_{res}_wholegenome.pdf"
	script:
		"scripts/plot_whole_genome.py"

# Figure 5


rule plot_microc_regions:
	input:
		"data/micro-c/merged/control/hic/control_{res}.hic"
	output: 
		expand("figures/figure_5_panels/{{res}}/{name}.pdf", 
			name=["twi", "sna", "if", "NetA", "sog", "Doc1", "pnr", "C15"])
	shell:
		"python scripts/plot_micro-c_regions.py {wildcards.res} "

rule plot_ep_loop_aggregates:
    input:
        loops = "data/loops/enhancer_promoter_loops/{loop_set}_loops.bedpe",
        hic = "data/hic/merged/{sample}/hic/{sample}_{res}.hic"
    output:
        plot = "figures/plot_ep_loop_aggregates/{loop_set}_{sample}_{res}_aggregate.pdf",
        matrix = "figures/plot_ep_loop_aggregates/{loop_set}_{sample}_{res}_aggregate.txt",
        aggmat = "figures/plot_ep_loop_aggregates/{loop_set}_{sample}_{res}_aggregate.hdf5"
    shell:
        "fanc aggregate --loops -tmp "
        "--colormap bwr --vmin -1 --vmax 1 --pixels 30 "
        "-m {output.matrix} -p {output.plot} "
        "{input.hic} {input.loops} {output.aggmat}"

rule write_ep_loop_strengths:
    input:
        loops = "data/loops/enhancer_promoter_loops/{loop_set}_loops.bedpe",
        hic = "data/hic/merged/{sample}/hic/{sample}_{res}.hic"
    output:
        "figures/plot_ep_loop_aggregates/{loop_set}_{sample}_{res}_loop_strengths.bedpe"
    shell:
        "python scripts/write_loop_strengths.py {input.loops} {input.hic} {output}"


rule plot_ep_loop_aggregates_microc:
    input:
        loops = "data/loops/enhancer_promoter_loops/{loop_set}_loops.bedpe",
        hic = "data/micro-c/merged/{sample}/hic/{sample}_{res}.hic"
    output:
        plot = "figures/plot_ep_loop_aggregates/{loop_set}_{sample}-micro-c_{res}_aggregate.pdf",
        matrix = "figures/plot_ep_loop_aggregates/{loop_set}_{sample}-micro-c_{res}_aggregate.txt",
        aggmat = "figures/plot_ep_loop_aggregates/{loop_set}_{sample}-micro-c_{res}_aggregate.hdf5"
    wildcard_constraints:
        sample = "gd7|control"
    shell:
        "fanc aggregate --loops -tmp "
        "--colormap bwr --vmin -1 --vmax 1 --pixels 30 "
        "-m {output.matrix} -p {output.plot} "
        "{input.hic} {input.loops} {output.aggmat}"

rule write_ep_loop_strengths_microc:
    input:
        loops = "data/loops/enhancer_promoter_loops/{loop_set}_loops.bedpe",
        hic = "data/micro-c/merged/{sample}/hic/{sample}_{res}.hic"
    output:
        "figures/plot_ep_loop_aggregates/{loop_set}_{sample}-micro-c_{res}_loop_strengths.bedpe"
    wildcard_constraints:
        sample = "gd7|control"
    shell:
        "python scripts/write_loop_strengths.py {input.loops} {input.hic} {output}"

rule plot_schematic:
    input:
        hic = "data/micro-c/merged/control/hic/control_500bp.hic",
        enh = "data/supplementary_tables/gd7_candidate_enhancers.bed"
    output:
        "figures/figure_5_panels/ep_loop_schematic.pdf"
    shell:
        """
        fancplot 3L:9,000,000-9,050,000 -o {output} --width 2 \
        --plot square {input.hic} --hide-x \
        --log -vmin 1e-3 -vmax 1e-1 --no-colorbar \
        --plot layer {input.enh} --hide-x \
        --plot gene external_data/flybase/dmel-all-r6.30.gtf.gz --hide-x \
        --no-labels --group-by gene_id --squash --aspect-ratio 0.1
        """

rule plot_microc_overview_figure:
    input: 
        "data/micro-c/merged/{merged_sample_name}/hic/{merged_sample_name}_10kb.hic",
        "data/micro-c/merged/{merged_sample_name}/hic/{merged_sample_name}_1kb.hic",
        "data/micro-c/merged/{merged_sample_name}/hic/{merged_sample_name}_50kb_masked.cor"
    output:
        "figures/figure_5_panels/{merged_sample_name}_overview_figure.pdf"
    shell:
        "python scripts/plot_micro-c_overview_figure.py {wildcards.merged_sample_name}"

## stg10 figure

rule plot_stg10_regions:
    input:
        "data/hic/merged/control-stg10/hic/control-stg10_2kb.hic"
    output: 
        expand("figures/figure_stg10_panels/{name}.pdf", 
            name=["twi", "sna", "if", "NetA", "sog", "Doc1", "pnr", "C15"])
    shell:
        "python scripts/plot_stg10_regions.py"
