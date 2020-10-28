# ###################
# # COMPARTMENTS
# ###################

rule custom_eigenvectors:
	input: 
		"data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}.hic"
	output:
		"data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}_eig1.bed",
 		"data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}_eig2.bed",
 		"data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}.cor"
	shell: "python scripts/eigenvectors.py {input}"

rule ref_eigenvector:
	input:
		"data/hic/merged/3-4h/hic/3-4h_{res}.cor"
	output:
		"data/hic/merged/3-4h/hic/3-4h_{res}_fanc_eigenvector.bed"
	shell: "fanc compartments -tmp -g {config[genome_fasta]} -v {output} {input}"

rule compartment_aggregates:
	input:
		hic = "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_50kb.hic",
		ab = "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_50kb.cor",
		ev = "data/hic/merged/3-4h/hic/3-4h_50kb_fanc_eigenvector.bed"
	output:
		plot = "figures/compartment_aggregates/{merged_sample_name}_50kb_comp_aggregate_{n_percentiles}.pdf",
		mat = "figures/compartment_aggregates/{merged_sample_name}_50kb_comp_aggregate_{n_percentiles}.npy.txt"
	params:
		percentiles = lambda wildcards: " ".join([str(x) for x in range(0, 101, 100 // int(wildcards.n_percentiles)) if x > 0])
	shell:
		"fanc compartments -tmp "
		"-p {params.percentiles} "
		"-v {input.ev} -e {output.plot} -m {output.mat} "
		"{input.hic} {input.ab}"


rule mask_heterochromatin:
	input: 
		hic = "data/{type}/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}.hic",
		mask = "external_data/modencode_white/pericentromeric_h3k9me3.bed"
	output: 
		"data/{type}/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}_masked.hic"
	script:
		"scripts/mask_heterochromatin.py"

rule custom_eigenvectors_masked:
	input: 
		"data/{type}/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}_masked.hic"
	output:
		"data/{type}/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}_masked_eig1.bed",
 		"data/{type}/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}_masked_eig2.bed",
 		"data/{type}/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}_masked.cor"
	shell: "python scripts/eigenvectors.py {input}"

rule ref_eigenvector_masked:
	input:
		"data/hic/merged/3-4h/hic/3-4h_{res}_masked.cor"
	output:
		"data/hic/merged/3-4h/hic/3-4h_{res}_masked_fanc_eigenvector.bed"
	shell: "fanc compartments -tmp -g {config[genome_fasta]} -v {output} {input}"

rule compartment_aggregates_masked:
	input:
		hic = "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_50kb_masked.hic",
		ab = "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_50kb_masked.cor",
		ev = "data/hic/merged/3-4h/hic/3-4h_50kb_masked_fanc_eigenvector.bed",
		ev2 = "data/compartments_by_gene_density/3-4h_50kb_masked_corrected_eigenvector.bed"
	output:
		plot = "figures/compartment_aggregates/{merged_sample_name}_50kb_masked_comp_aggregate_{n_percentiles}.pdf",
		mat = "figures/compartment_aggregates/{merged_sample_name}_50kb_masked_comp_aggregate_{n_percentiles}.npy.txt"
	params:
		percentiles = lambda wildcards: " ".join([str(x) for x in range(0, 101, 100 // int(wildcards.n_percentiles)) if x > 0])
	shell:
		"fanc compartments -tmp "
		"-p {params.percentiles} "
		"-v {input.ev2} -e {output.plot} -m {output.mat} "
		"{input.hic} {input.ab}"


rule compartment_aggregates_masked_collapsed:
	input:
		hic = "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_50kb_masked.hic",
		ab = "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_50kb_masked.cor",
		ev = "data/hic/merged/3-4h/hic/3-4h_50kb_masked_fanc_eigenvector.bed",
		ev2 = "data/compartments_by_gene_density/3-4h_50kb_masked_corrected_eigenvector.bed"
	output:
		plot = "figures/compartment_aggregates/{merged_sample_name}_50kb_masked_collapsed_comp_aggregate_{n_percentiles}.pdf"
	params:
		n = lambda wildcards: int(wildcards.n_percentiles)
	script:
		"scripts/plot_compartment_aggregates_collapsed.py"

