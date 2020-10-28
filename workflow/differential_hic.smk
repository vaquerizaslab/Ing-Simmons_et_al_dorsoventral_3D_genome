rule export_uncorrected_hic_as_text_for_comparisons:
    input:
        ancient("data/hic/{sample}/hic/{sample}_{res}.hic")
    output:
        regions = "data/hic/{sample}/hic/{sample}_{res}_uncorrected_regions.bed",
        matrix = "data/hic/{sample}/hic/{sample}_{res}_uncorrected_matrix.txt"
    shell:
        "./fan-c/bin/fanc dump --uncorrected {input} {output.matrix} {output.regions}"


rule export_bias_for_ACCOST:
    input:
        ancient("data/hic/{sample}/hic/{sample}_{res}.hic")
    output:
        "data/hic/{sample}/hic/{sample}_{res}_bias.txt"
    shell:
        "python3 ./scripts/export_bias.py {input} {output}"

rule diffHic:
    input:
        "data/hic/{sample1}_Rep1/hic/{sample1}_Rep1_{res}_uncorrected_regions.bed",
        "data/hic/{sample1}_Rep1/hic/{sample1}_Rep1_{res}_uncorrected_matrix.txt",
        "data/hic/{sample2}_Rep1/hic/{sample2}_Rep1_{res}_uncorrected_regions.bed",
        "data/hic/{sample2}_Rep1/hic/{sample2}_Rep1_{res}_uncorrected_matrix.txt",
        "data/hic/{sample1}_Rep2/hic/{sample1}_Rep2_{res}_uncorrected_regions.bed",
        "data/hic/{sample1}_Rep2/hic/{sample1}_Rep2_{res}_uncorrected_matrix.txt",
        "data/hic/{sample2}_Rep2/hic/{sample2}_Rep2_{res}_uncorrected_regions.bed",
        "data/hic/{sample2}_Rep2/hic/{sample2}_Rep2_{res}_uncorrected_matrix.txt"
    output:
        report = "scripts/diffHic_{sample1}_vs_{sample2}_{res}.html",
        results = "data/diffHic/{sample1}_vs_{sample2}_{res}_diffHic_results.txt"
    shell:
        """
        R -e \"rmarkdown::render('/home/research/vaquerizas/liz/dorsal_ventral/fan-c/scripts/diffhic.Rmd', \
            output_file='/home/research/vaquerizas/liz/dorsal_ventral/fan-c/{output.report}',\
            params = list(sample1 = '{wildcards.sample1}', sample2 = '{wildcards.sample2}', res = '{wildcards.res}'),\
            knit_root_dir='/home/research/vaquerizas/liz/dorsal_ventral/fan-c/scripts/')\"
        """
        

rule diffHic_analysis:
    input:
        expand("data/diffHic/{sample1}_vs_{sample2}_{res}_diffHic_results.txt",
            sample1 = ["control-nc14"], res = ["5kb", "10kb"],
            sample2 = ["gd7-nc14", "Tollrm910-nc14", "Toll10B-nc14"])
    output:
        expand("data/diffHic/{sample1}_vs_{sample2}_{res}_top_regions.bed",
            sample1 = ["control-nc14"], res = ["5kb", "10kb"],
            sample2 = ["gd7-nc14", "Tollrm910-nc14", "Toll10B-nc14"]),
        expand("data/diffHic/{sample1}_vs_{sample2}_{res}_{direction}_interactions.bedpe",
            sample1 = ["control-nc14"], res = ["5kb", "10kb"],
            sample2 = ["gd7-nc14", "Tollrm910-nc14", "Toll10B-nc14"], direction = ["up", "down"]),
        report = "scripts/diffhic_analysis.html"
    shell:
        """
        R -e \"rmarkdown::render('/home/research/vaquerizas/liz/dorsal_ventral/fan-c/scripts/diffhic_analysis.Rmd', \
            output_file='/home/research/vaquerizas/liz/dorsal_ventral/fan-c/{output.report}',\
            knit_root_dir='/home/research/vaquerizas/liz/dorsal_ventral/fan-c/scripts/')\"
        """

rule plot_diffhic_regions:
    input:
        "data/diffHic/{sample1}_vs_{sample2}_{res}_top_regions_with_anchors.bed"
    output:
        "figures/diffHic/{sample1}_vs_{sample2}_{res}_top_regions.pdf"
    shell:
        "python scripts/plot_diffhic_regions.py {input} {output}"

rule plot_diffhic_aggregates:
    input:
        hic = "data/hic/merged/{sample}/hic/{sample}_{res}.hic",
        loops = "data/diffHic/{sample1}_vs_{sample2}_{res}_{direction}_all_interactions.bedpe"
    output:
        plot = "figures/diffHic/{sample1}_vs_{sample2}_{res}_{direction}_{sample}_aggregate.pdf",
        matrix = "figures/diffHic/{sample1}_vs_{sample2}_{res}_{direction}_{sample}_aggregate.txt",
        aggmat = "figures/diffHic/{sample1}_vs_{sample2}_{res}_{direction}_{sample}_aggregate.hdf5"
    shell:
        "./fan-c/bin/fanc aggregate --loops -tmp "
        "--colormap bwr --vmin -1 --vmax 1 --pixels 100 "
        "-m {output.matrix} -p {output.plot} "
        "{input.hic} {input.loops} {output.aggmat}"


        
