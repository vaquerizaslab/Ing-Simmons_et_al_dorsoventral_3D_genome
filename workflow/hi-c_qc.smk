rule cis_trans_ratio:
    input:
        "data/{type}/{sample_name}/hic/{sample_name}.hic"
    output:
        "data/{type}/{sample_name}/{sample_name}_cis_trans_ratio.txt"
    shell: "fanc cis_trans -o {output} {input}"

rule cis_local_ratio:
    input:
        "data/{type}/{sample_name}/hic/{sample_name}.hic"
    output:
        "data/{type}/{sample_name}/{sample_name}_cis_local_ratio.txt"
    params:
        threshold = 20000
    shell: "python scripts/get_cis_local_ratio.py {input} {params.threshold} {output}"

rule get_hic_stats:
    input: 
        "data/{type}/{sample_name}/hic/{sample_name}.hic"
    output:
        "data/{type}/{sample_name}/{sample_name}_hic_stats.txt"
    shell: "python scripts/get_stats.py data/{wildcards.type}/{wildcards.sample_name}/pairs/ {output}"

def get_fastq(wildcards):
    fastq_files = metadata.iloc[:,3:].values.ravel()
    fastq_files = fastq_files[~pd.isnull(fastq_files)].tolist()
    fastq_files = [os.path.join("data/fastq", f) for f in fastq_files]

    microc_fastq_files = microc_metadata.iloc[:,3:].values.ravel()
    microc_fastq_files = microc_fastq_files[~pd.isnull(microc_fastq_files)].tolist()
    microc_fastq_files = [os.path.join("data/fastq", f) for f in microc_fastq_files]

    return fastq_files + microc_fastq_files

def get_bam(wildcards):
    def get_bam_names(row):
        fastqs = row.iloc[3:]
        fastqs = fastqs[~pd.isnull(fastqs)].tolist()
        x = [os.path.join("data/hic", row.condition + "_" + row.replicate,  "sam", re.sub(".fastq.gz|.fq.gz", ".bam", x)) for x in fastqs]
        return(x)

    bam_files = [get_bam_names(row) for index, row in metadata.iterrows() ]
    bam_files = [item for sublist in bam_files for item in sublist]

    microc_bam_files = [get_bam_names(row) for index, row in microc_metadata.iterrows() ]
    microc_bam_files = [item for sublist in microc_bam_files for item in sublist]    

    return bam_files + microc_bam_files

rule get_alignment_stats_bwa:
    input:
        unpack(get_fastq),
        unpack(get_bam)
    output:
        "data/sample_alignment_stats_bwa.txt"
    threads: 16 
    script: "scripts/alignment_stats_bwa.R"

rule plot_stats:
    input:
        "data/sample_alignment_stats_bwa.txt",
        expand("data/hic/{sample}/{sample}_hic_stats.txt", sample = SAMPLE_NAMES),
        expand("data/hic/{sample}/{sample}_cis_trans_ratio.txt", sample = SAMPLE_NAMES),
        expand("data/hic/{sample}/{sample}_cis_local_ratio.txt", sample = SAMPLE_NAMES),
        expand("data/micro-c/{sample}/{sample}_hic_stats.txt", sample = microc_metadata.index),
        expand("data/micro-c/{sample}/{sample}_cis_trans_ratio.txt", sample = microc_metadata.index),
        expand("data/micro-c/{sample}/{sample}_cis_local_ratio.txt", sample = microc_metadata.index)
    output:
        "scripts/plot_stats.html"
    script: "scripts/plot_stats.Rmd"


