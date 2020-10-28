fastq_dict = {s: ["data/fastq/{}".format(file) for file in metadata.loc[s, metadata.columns.str.startswith('run')] if not pd.isnull(file)] for s in metadata.index}
bam_dict = {}

for sample, fastqs in fastq_dict.items():
    bam_dict[sample] = [os.path.join("data/hic", sample,  "sam", re.sub(".fastq.gz|.fq.gz", ".bam", os.path.basename(x))) for x in fastqs]

pairs_dict = {}
for sample, fastqs in fastq_dict.items():
    n_runs = len(fastqs)/2
    if not n_runs.is_integer():
        raise AssertionError("Odd number of fastq files provided for {}!", sample)
    else:
        pairs_dict[sample] = [os.path.join("data/hic", sample, "pairs", sample + "_" + str(i) + ".pairs") for i in range(int(n_runs))]

for sample in metadata.index:
    rule:
        input:
            fastq_dict[sample]
        output:
            bam_dict[sample]
        threads: 16
        params:
            name = sample,
            output_folder = f"data/hic/{sample}/sam/"
        shell:
            "fanc map -t {threads} -tmp "
            " -q 3 "
            "{input} {config[genome_idx]} {params.output_folder} "

rule make_fragments:
    output:
        "data/dm6_RE_fragments.bed"
    params:
        chrs = ",".join(chromosomes),
        genome = config["genome_fasta"],
        enzyme = config["enzyme"]
    shell:
        "fanc fragments -c {params.chrs} {params.genome} {params.enzyme} {output}"

for sample in metadata.index:
    rule:
        input:
            bams = bam_dict[sample],
            genome = "data/dm6_RE_fragments.bed"
        output:
            pairs_dict[sample]
        threads: 4
        params:
            name = sample,
            output_folder = f"data/hic/{sample}"
        shell:
            "fanc auto -t {threads} -tmp "
            "-g {input.genome} "
            " -q 3 "
            "--le-inward-cutoff 1000 --le-outward-cutoff 1000 "
            "--no-hic -n {params.name} "
            "{input.bams} {params.output_folder} "

for sample in metadata.index:
    rule:
        input:
            pairs_dict[sample]
        output:
            f"data/hic/{sample}/hic/{sample}.hic"
        shell:
            "fanc hic -tmp {input} {output}"

rule process_hic:
    input:
        "data/hic/{sample}/hic/{sample}.hic"
    output:
        out = "data/hic/{sample}/hic/{sample}_{res}.hic",
        stats = "data/hic/{sample}/hic/{sample}_{res}_filter_stats.txt",
        stats_plot = "data/hic/{sample}/hic/{sample}_{res}_filter_stats.pdf",
    threads: 8
    shell:
        "fanc hic -t {threads} -tmp "
        "-b {wildcards.res} --low-coverage-auto --kr-correct "
        "-s {output.stats} --statistics-plot {output.stats_plot} "
        "{input} {output.out}"


rule export_for_pca:
    input:
        "data/hic/{sample}/hic/{sample}_100kb.hic"
    output:
        expand("data/hic/{{sample}}/hic/{{sample}}_100kb_{chr}.npy.txt",
            chr = ["2L", "2R", "3L", "3R", "X", "all"])
    params:
        upper = "50",
        lower = "1"
    shell: "python scripts/export_npy_matrix.py {input} {params.upper} {params.lower}"

rule fanc_pca:
    input: 
        expand("data/hic/{sample}/hic/{sample}_100kb.hic", sample = metadata.index)
    output:
        pca_file = "data/hic/fanc_pca_100kb.pca",
        plot = "data/hic/fanc_pca_100kb.pdf"
    shell:
        "fanc pca -p {output.plot} "
        "--min-distance 100kb --max-distance 1Mb "
        "{input} {output.pca_file} "
