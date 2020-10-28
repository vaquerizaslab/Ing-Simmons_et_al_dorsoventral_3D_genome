###################
# CHESS
###################

rule export_hic_as_text:
  input:
    ancient("data/hic/merged/{sample}/hic/{sample}_{res}.hic")
  output:
    regions = "data/hic/merged/{sample}/hic/{sample}_{res}_regions.bed",
    matrix = "data/hic/merged/{sample}/hic/{sample}_{res}_matrix.txt"
  shell:
    "fanc dump {input} {output.matrix} {output.regions}"

rule chess_make_pairs:
  output: 
    "data/chess/dm6_pairs_{window}x_{res}.bedpe"
  params:
    window_bp = lambda wildcards: 1000 * int(wildcards.window) * int(re.sub('kb', '', wildcards.res)),
    step_bp = lambda wildcards: 1000 * int(re.sub('kb', '', wildcards.res))
  shell: 
      "chess pairs ../dm6_chrom_sizes_sanitized.txt {params.window_bp} {params.step_bp} {output} --file-input "

rule chess_sim:
  input:
    matrix1 = "data/hic/merged/{sample1}/hic/{sample1}_{res}.hic",
    matrix2 = "data/hic/merged/{sample2}/hic/{sample2}_{res}.hic",
    regions1 = "data/hic/merged/{sample1}/hic/{sample1}_{res}_regions.bed",
    regions2 = "data/hic/merged/{sample2}/hic/{sample2}_{res}_regions.bed",
    pairs = "data/chess/dm6_pairs_{window}x_{res}.bedpe"
  output:
    "data/chess/{sample1}_vs_{sample2}/genome_scan_{window}x_{res}.txt"
  threads: 8
  shell:
    "chess sim "
    # "--reference-regions {input.regions1} --query-regions {input.regions2} "
    "{input.matrix1} {input.matrix2}  "
    "{input.pairs} {output} "
    "-p {threads}"

rule analyse_chess:
  input:
    expand(expand("data/chess/{sample1}_vs_{sample2}/genome_scan_{{window}}x_{{res}}.txt",
      zip, sample1 = REF_SAMPLES, sample2 = QUERY_SAMPLES),
      res = ["5kb", "10kb", "25kb", "50kb"], window = ["100", "150"])
  output:
    "scripts/analyse_chess_comparisons_merged.html",
    # N.B. Rmarkdown rendering with script only works with a single output file
    expand("data/chess/{res}_{window}x_signif_regions.bed",
      res = ["5kb", "10kb", "25kb", "50kb"], window = ["100", "150"]),
    expand(expand("data/chess/{sample1}_vs_{sample2}/{{res}}_{{window}}x_signif_regions.bedpe",
      zip, sample1 = REF_SAMPLES, sample2 = QUERY_SAMPLES),
      res = ["5kb", "10kb", "25kb", "50kb"], window = ["100", "150"])
  # script:
  #   "scripts/analyse_chess_comparisons_merged.Rmd"
  shell:
    "R -e \"rmarkdown::render('/home/research/vaquerizas/liz/dorsal_ventral/fan-c/scripts/analyse_chess_comparisons_merged.Rmd', "
    "output_file='/home/research/vaquerizas/liz/dorsal_ventral/fan-c/scripts/analyse_chess_comparisons_merged.html', "
    "knit_root_dir='/home/research/vaquerizas/liz/dorsal_ventral/fan-c/scripts/')\""

# rule chess_extract:
#   input:
#     regions1 = "data/hic/merged/{sample1}/hic/{sample1}_{res}_regions.bed",
#     matrix1 = "data/hic/merged/{sample1}/hic/{sample1}_{res}_matrix.txt",
#     regions2 = "data/hic/merged/{sample2}/hic/{sample2}_{res}_regions.bed",
#     matrix2 = "data/hic/merged/{sample2}/hic/{sample2}_{res}_matrix.txt",
#     pairs = "data/chess/{sample1}_vs_{sample2}/{res}_{window}_signif_regions.bedpe"
#   output:
#     expand("data/chess/{{sample1}}_vs_{{sample2}}/{{res}}_{{window}}_features/{direction}_features.tsv",
#       direction = ['gained', 'lost'])
#   params:
#     outdir = "data/chess/{sample1}_vs_{sample2}/{res}_{window}_features/"
#   shell:
#     "chess extract "
#     "{input.pairs} "
#     "{input.matrix1} {input.regions1} {input.matrix2} {input.regions2} "
#     "{params.outdir} "

# rule chess_crosscorrelate:
#   input:
#     pairs = "data/chess/{sample1}_vs_{sample2}/{res}_{window}x_signif_regions.bedpe",
#     expand("data/chess/{sample1}_vs_{sample2}/{res}_{window}_features/{direction}_features.tsv",
#       direction = ['gained', 'lost'])
#   output:
    
#   shell:
#     "chess crosscorrelate "
#     "{input.pairs} "
#     "{input.matrix1} {input.regions1} {input.matrix2} {input.regions2} "
#     "{output} "

rule plot_chess_regions:
  input:
    regions ="data/chess/{res}_{window}x_signif_regions.bed",
    fc = expand("data/hic/merged/{sample}/hic/foldchange_control-nc14_{sample}_{res}.hic",
      sample = ["gd7-nc14", "Tollrm910-nc14", "Toll10B-nc14"], res = ["5kb", "10kb", "25kb", "50kb"])
  output:
    "figures/chess/{res}_{window}x_signif_regions.pdf"
  shell:
    "python scripts/plot_chess_regions2.py {input.regions} {output}"
