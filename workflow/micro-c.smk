microc_samples = ["gd7", "control"]
microc_resolutions = ["500kb", "100kb", "50kb", "25kb", "10kb", "5kb", "2kb", "1kb", "500bp", "100bp"]
microc_possible_replicates = ["Rep1", "Rep2", "Rep3", "Rep4"]


microc_fastq_files = [re.sub(".gz", "", f) for f in
               list(itertools.chain(*microc_metadata.loc[:, microc_metadata.columns.str.startswith('run')].values)) if not pd.isnull(f) ]

fastq_dict = {s: ["data/fastq/{}".format(file) for file in microc_metadata.loc[s, microc_metadata.columns.str.startswith('run')] if not pd.isnull(file)] for s in microc_metadata.index}

bam_dict = {}

for sample, fastqs in fastq_dict.items():
    bam_dict[sample] = [os.path.join("data/micro-c", sample,  "sam", re.sub(".fastq.gz|.fq.gz", ".bam", os.path.basename(x))) for x in fastqs]

bam_dict2 = {}

for sample, fastqs in fastq_dict.items():
    bam_dict2[sample] = [os.path.join("data/micro-c", sample,  "sam", re.sub(".fastq.gz|.fq.gz", "_sorted.bam", os.path.basename(x))) for x in fastqs]

pairs_dict = {}
for sample, fastqs in fastq_dict.items():
    n_runs = len(fastqs)/2
    if not n_runs.is_integer():
        raise AssertionError("Odd number of fastq files provided for {}!", sample)
    else:
        pairs_dict[sample] = [os.path.join("data/micro-c", sample, "pairs", sample + "_" + str(i) + ".pairs") for i in range(int(n_runs))]

for sample in microc_metadata.index:
    rule:
        input:
            fastq_dict[sample]
        output:
            bam_dict[sample]
        threads: 16
        params:
            name = sample,
            output_folder = f"data/micro-c/{sample}/sam/"
        shell:
            "fanc map -t {threads} -tmp "
            " -q 3 "
            "{input} {config[genome_idx]} {params.output_folder} "

rule get_microc_alignment_stats:
    input:
        list(fastq_dict.values()),
        list(bam_dict.values())
    output:
        "data/micro-c/sample_alignment_stats_bwa.txt"
    threads: 16 
    script: "scripts/micro-c_alignment_stats.R"

rule bin_genome:
    output:
        "data/dm6_100bp_fragments.bed"
    params:
        chrs = ",".join(chromosomes),
        genome = config["genome_fasta"]
    shell:
        "fanc fragments -c {params.chrs} {params.genome} 100 {output}"

rule sort_sam:
    input:
        "data/micro-c/{sample}/sam/{read}.bam"
    output:
        "data/micro-c/{sample}/sam/{read}_sorted.bam"
    threads: 4
    shell:
        "fanc sort-sam -t {threads} -tmp {input} {output}"

for sample in microc_metadata.index:
    rule:
        input:
            bams = bam_dict2[sample],
            genome = "data/dm6_100bp_fragments.bed"
        output:
            f"data/micro-c/{sample}/pairs/{sample}.pairs"
        shell:
            "fanc pairs -tmp "
            "-g {input.genome} "
            "--filter-unmappable --filter-multimapping --filter-quality 3 "
            "{input.bams} {output} "

rule filter_microc_pairs:
    input:
        pairs = ancient("data/micro-c/{sample}/pairs/{sample}.pairs"),
        genome = "data/dm6_100bp_fragments.bed"
    output:
        stats_plot = "data/micro-c/{sample}/pairs/{sample}_filtered_stats.png",
        redist_plot = "data/micro-c/{sample}/pairs/{sample}_filtered_re_dist.png",
        ligation_plot = "data/micro-c/{sample}/pairs/{sample}_filtered_ligation_error.png",
    shell:
        "fanc pairs -tmp "
        "-g {input.genome} "
        "--filter-inward 50 --filter-outward 50 "
        "--filter-self-ligations --filter-pcr-duplicates 1 "
        "--statistics-plot {output.stats_plot} "
        "--re-dist-plot {output.redist_plot} "
        "--ligation-error-plot {output.ligation_plot} "
        "{input.pairs} "

rule make_microc:
    input:
        pairs = "data/micro-c/{sample}/pairs/{sample}.pairs",
        ligation_plot = "data/micro-c/{sample}/pairs/{sample}_filtered_ligation_error.png",
    output:
        "data/micro-c/{sample}/hic/{sample}.hic"
    shell:
        "fanc hic -tmp {input.pairs} {output}"

def check_microc_reps(name):
    possible_names = snakemake.io.expand("{name}_{rep}", name = [name], rep = microc_possible_replicates)
    possible_paths = [f"data/micro-c/{name}/hic/{name}.hic" for name in possible_names]
    paths = [path for path in possible_paths if os.path.exists(path)]
    return(paths)

microc_dict = {name : check_microc_reps(name) for name in microc_samples}

for sample in microc_samples:
    rule:
        input:
            microc_dict[sample]
        output:
            f"data/micro-c/merged/{sample}/hic/{sample}.hic"
        threads: 1
        shell:
            "fanc hic -tmp {input} {output}"    

rule fanc_process_microc_merged:
    input: 
        "data/micro-c/merged/{sample}/hic/{sample}.hic"
    output:
        out = "data/micro-c/merged/{sample}/hic/{sample}_{res}.hic",
        stats = "data/micro-c/merged/{sample}/stats/{sample}_{res}_filter_stats.txt",
        stats_plot = "data/micro-c/merged/{sample}/hic/{sample}_{res}_filter_stats.pdf"
    threads: 8
    shell:
       "fanc hic -t {threads} -tmp "
        "-b {wildcards.res} --low-coverage-auto "
        "--normalise --norm-method ice "
        "-s {output.stats} --statistics-plot {output.stats_plot} "
        "{input} {output.out}"

rule make_cooler_microc:
    input:
        "data/micro-c/merged/{merged_sample_name}/hic/{merged_sample_name}_{res}.hic"
    output:
        "data/micro-c/cooler_files/{merged_sample_name}_{res}.mcool"
    threads: 1
    shell:
        "fanc to-cooler {input} {output} "

rule export_marginals_microc:
    input:
        "data/micro-c/merged/{sample}/hic/{sample}_{res}.hic"
    output:
        "data/micro-c/merged/{sample}/hic/{sample}_{res}_marginals.bed"
    shell: 
        "python scripts/export_marginals.py {input} {output}"

# rule calculate_resolution:
#     input:
#         expand("data/micro-c/merged/{sample}/hic/{sample}_{res}_marginals.bed",
#             sample = microc_samples, res=microc_resolutions)
#     output:
#         "scripts/analyse_micro-c_resolution.html"
#     script:
#         "scripts/analyse_micro-c_resolution.Rmd"
