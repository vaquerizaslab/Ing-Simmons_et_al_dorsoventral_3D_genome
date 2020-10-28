# ###################
# # INSULATION
# ###################

rule insulation_score:
	input: 
		"data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}.hic"
	output:
		"data/boundaries/{merged_sample_name}_{res}.ii",
		expand("data/boundaries/{{merged_sample_name}}_{{res}}_{w}.bw",
			w = [4, 6, 8, 10])
	shell: "python scripts/insulation_score.py {input}"

rule call_boundaries:
	input: 
		"data/boundaries/{merged_sample_name}_{res}.ii"
	output:
		"data/boundaries/{merged_sample_name}_{res}_{w}bins_d3_boundaries.bed"
	params:
		window = lambda wildcards: int(re.findall('([0-9]+)kb', wildcards.res)[0]) * int(wildcards.w) * 1000
	shell:
		"fanc boundaries -w {params.window} -d 3 {input} {output}"

rule filter_boundaries:
	input:
		"data/boundaries/{merged_sample_name}_{res}_{w}bins_d3_boundaries.bed"
	output:
		"data/boundaries/{merged_sample_name}_{res}_{w}bins_d3_boundaries_filtered0.7.bed"
	params:
		threshold = 0.7
	shell:
		"Rscript scripts/filter_boundaries.R {input} {params.threshold}"

rule make_boundaries_and_domains:
	input:
		expand("data/boundaries/{merged_sample_name}_{res}_{w}bins_d3_boundaries_filtered0.7.bed",
			merged_sample_name = merged_sample_names, res = ["2kb", "5kb"], w = [4, 6, 8, 10])
	output:
		expand("data/boundaries/{merged_sample_name}_final_boundaries.bed", merged_sample_name = merged_sample_names),
		expand("data/boundaries/{merged_sample_name}_paired_boundary_domains.bed", merged_sample_name = merged_sample_names),
		"scripts/domain_boundaries.html"
	shell:
		"R -e \"rmarkdown::render('/home/research/vaquerizas/liz/dorsal_ventral/fan-c/scripts/domain_boundaries.Rmd', "
		"output_file='/home/research/vaquerizas/liz/dorsal_ventral/fan-c/scripts/domain_boundaries.html', "
		"knit_root_dir='/home/research/vaquerizas/liz/dorsal_ventral/fan-c/scripts/')\""

rule make_boundaries_hitile:
	input:
		"data/boundaries/{merged_sample_name}_final_boundaries.bed"
	output:
		"data/boundaries/{merged_sample_name}_final_boundaries.hitile"
	shell:
		"clodius aggregate bedfile "
		" --no-header "
		" --chromsizes-filename ../dm6_chrom_sizes_sanitized.txt "
		" -o {output} {input} "

rule make_domains_hitile:
	input:
		"data/boundaries/{merged_sample_name}_paired_boundary_domains.bed"
	output:
		"data/boundaries/{merged_sample_name}_paired_boundary_domains.hitile"
	shell:
		"clodius aggregate bedfile "
		" --no-header "
		" --chromsizes-filename ../dm6_chrom_sizes_sanitized.txt "
		" -o {output} {input} "

rule plot_domain_aggregates:
	input:
		domains = "data/boundaries/3-4h_paired_boundary_domains.bed",
		hic = "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}.hic"
	output:
		plot = "figures/plot_domain_aggregates/3-4h_paired_domains_{merged_sample_name}_{res}_aggregate.pdf",
		matrix = "figures/plot_domain_aggregates/3-4h_paired_domains_{merged_sample_name}_{res}_aggregate.txt",
		aggmat = "figures/plot_domain_aggregates/3-4h_paired_domains_{merged_sample_name}_{res}_aggregate.hdf5"
	shell:
		"fanc aggregate --tads -tmp "
		"--colormap bwr --vmin -1 --vmax 1 "
		"-m {output.matrix} -p {output.plot} "
		"{input.hic} {input.domains} {output.aggmat}"

rule plot_boundary_aggregates:
	input:
		boundaries = "data/boundaries/3-4h_final_boundaries.bed",
		hic = ancient("data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}.hic")
	output:
		plot = "figures/plot_boundary_aggregates/3-4h_boundaries_{merged_sample_name}_{res}_aggregate.pdf",
		matrix = "figures/plot_boundary_aggregates/3-4h_boundaries_{merged_sample_name}_{res}_aggregate.txt",
		aggmat = "figures/plot_boundary_aggregates/3-4h_boundaries_{merged_sample_name}_{res}_aggregate.hdf5"
	shell:
		"fanc aggregate -tmp "
		"--window 100kb -e --log "
		"--colormap bwr --vmin -1 --vmax 1 "
		"-m {output.matrix} -p {output.plot} "
		"{input.hic} {input.boundaries} {output.aggmat}"
