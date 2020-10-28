
rule process_published_loop_sets:
	output: 
		expand("../external_data/{publication}/{publication}_loops.bedpe", 
			publication = ["cubenas-potts_2016", "eagen_2017", "stadler_2017"])
	script: 
		"../external_data/process_loops.R"

rule plot_published_loop_aggregates:
	input:
		loops = "../external_data/{publication}/{publication}_loops.bedpe",
		hic = ancient("data/hic/merged/{sample}/hic/{sample}_{res}.hic")
	output:
		plot = "figures/plot_published_loop_aggregates/{publication}_{sample}_{res}_aggregate.pdf",
		matrix = "figures/plot_published_loop_aggregates/{publication}_{sample}_{res}_aggregate.txt",
		aggmat = "figures/plot_published_loop_aggregates/{publication}_{sample}_{res}_aggregate.hdf5"
	shell:
		"fanc aggregate --loops -tmp "
		"--colormap bwr --vmin -1 --vmax 1 --pixels 30 "
		"-m {output.matrix} -p {output.plot} "
		"{input.hic} {input.loops} {output.aggmat}"

rule write_loop_strengths:
	input:
		loops = "../external_data/{publication}/{publication}_loops.bedpe",
		hic = "data/hic/merged/{sample}/hic/{sample}_{res}.hic"
	output:
		"figures/plot_published_loop_aggregates/{publication}_{sample}_{res}_loop_strengths.bedpe"
	shell:
		"python scripts/write_loop_strengths.py {input.loops} {input.hic} {output}"
