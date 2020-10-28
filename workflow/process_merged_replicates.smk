###################
# MERGE HI-C FILES
###################

def check_reps(name):
    possible_names = snakemake.io.expand("{name}_{rep}", name = [name], rep = possible_replicates)
    possible_paths = [f"data/hic/{name}/hic/{name}.hic" for name in possible_names]
    paths = [path for path in possible_paths if os.path.exists(path)]
    return(paths)

hic_dict = {name : check_reps(name) for name in merged_sample_names}

for sample in merged_sample_names:
    rule:
        input:
            hic_dict[sample]
        output:
            f"data/hic/merged/{sample}/hic/{sample}.hic"
        threads: 1
        shell:
            "fanc hic -tmp {input} {output}"    

rule fanc_process_hic_merged:
    input: 
        "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}.hic"
    output:
        out = "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}.hic",
        stats = "data/hic/merged/{merged_sample_name}/stats/{merged_sample_name}_{res}_filter_stats.txt",
        stats_plot = "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}_filter_stats.pdf"
    threads: 8
    shell:
        "fanc hic -t {threads} -tmp "
        "-b {wildcards.res} --low-coverage-auto --kr-correct "
        "-s {output.stats} --statistics-plot {output.stats_plot} "
        "{input} {output.out}"

rule balance_whole_genome:
    input: 
        "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}.hic"
    output:
        "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}_wholegenomebalanced.hic"
    threads: 8
    shell:
        "fanc hic -t {threads} -tmp "
        "-b {wildcards.res} --low-coverage-auto --kr-correct --whole-matrix "
        "{input} {output}"



rule make_cooler:
    input:
        "data/hic/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}.hic"
    output:
        "data/hic/cooler_files/{merged_sample_name}_{res}.mcool"
    threads: 1
    shell:
        "fanc to-cooler {input} {output} "

rule make_diff_hic:
    input: 
        "data/hic/merged/{sample}/hic/{sample}_{res}.hic"
    output:
        "data/hic/merged/{sample}/hic/diff_control-nc14_{sample}_{res}.hic"
    shell: 
        "fanc compare -c difference -tmp --no-scale "
        "data/hic/merged/control-nc14/hic/control-nc14_{wildcards.res}.hic {input} {output}"



## Additional QC

rule calc_distance_decay:
    input: 
        "data/{type}/merged/{sample}/hic/{sample}_{res}.hic"
    output:
        "data/{type}/merged/{sample}/hic/{sample}_{res}_expected_values_all.txt",
        "data/{type}/merged/{sample}/hic/{sample}_{res}_expected_values_per_chrom.txt"
    shell: 
        "python scripts/calc_distance_decay.py {input} "
        "data/{wildcards.type}/merged/{wildcards.sample}/hic/{wildcards.sample}_{wildcards.res}"

rule analyse_distance_decay:
    input:
        expand("data/hic/merged/{sample}/hic/{sample}_{res}_expected_values_all.txt",
            sample = merged_sample_names, res = RESOLUTIONS),
        expand("data/hic/merged/{sample}/hic/{sample}_{res}_expected_values_per_chrom.txt",
            sample = merged_sample_names, res = RESOLUTIONS)
    output:
        "scripts/analyse_distance_decay.html"
    script:
        "analyse_distance_decay.Rmd"


rule export_marginals:
    input:
        "data/hic/merged/{sample}/hic/{sample}_{res}.hic"
    output:
        "data/hic/merged/{sample}/hic/{sample}_{res}_marginals.bed"
    shell: 
        "python scripts/export_marginals.py {input} {output}"

rule calculate_resolution:
    input:
        expand("data/hic/merged/{sample}/hic/{sample}_{res}_marginals.bed",
            sample = merged_sample_names, res=RESOLUTIONS)
    output:
        "scripts/analyse_resolution.html"
    script:
        "scripts/analyse_resolution.Rmd"

