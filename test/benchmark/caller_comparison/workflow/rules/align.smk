localrules: index_alignment, alignment_stats, pool_samples

def get_samples(wildcards):
    return config["samples"][wildcards.sample]

include: "align/pbmm2.smk"
include: "align/minimap2.smk"
include: "align/lra.smk"

rule index_alignment:
    input:
        "{name}.bam"
    output:
        "{name}.bam.bai"
    threads: 1
    resources:
        io_gb = 100
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule alignment_stats:
    input:
        bam = expand("pipeline/alignments/{sample}.{{aligner}}.bam", sample=config["samples"]),
        bai = expand("pipeline/alignments/{sample}.{{aligner}}.bam.bai", sample=config["samples"])
    output:
        "pipeline/alignment_stats/alignment_stats.{aligner}.txt"
    resources:
        io_gb = 100
    log:
        "pipeline/logs/alignment_stats/alignment_stats.{aligner}.log"
    shell:
        "python3 workflow/scripts/alignment_stats.py -o {output} {input.bam} 2> {log}"

rule pool_samples:
    input:
        expand("pipeline/alignments/{sample}.{{aligner}}.bam", sample=config["samples"])
    output:
        "pipeline/alignment_pooled/pooled.{aligner}.bam"
    resources:
        io_gb = 100
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools merge {output} {input}"

rule subsample_alignments_0:
    input:
        bam = "pipeline/alignment_pooled/pooled.{aligner}.bam"
    output:
        expand("pipeline/alignment_pooled/pooled.subsampled.{fraction}.{{aligner}}.bam", fraction=range(30, 100, 30))
    threads: 10
    resources:
        mem_mb = 400000,
        time_min = 1000,
        io_gb = 100
    params:
        tmpdir = "1000",
        outdir = "pipeline/alignment_pooled/"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bash workflow/scripts/subsample.sh {input.bam} 30 90 30 {threads} {params.outdir}"

rule subsample_alignments_1:
    input:
        bam = "pipeline/alignment_pooled/pooled.{aligner}.bam"
    output:
        expand("pipeline/alignment_pooled/pooled.subsampled.{fraction}.{{aligner}}.bam", fraction=range(10, 100, 30))
    threads: 10
    resources:
        mem_mb = 400000,
        time_min = 1000,
        io_gb = 100
    params:
        tmpdir = "1000",
        outdir = "pipeline/alignment_pooled/"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bash workflow/scripts/subsample.sh {input.bam} 10 90 30 {threads} {params.outdir}"

rule subsample_alignments_2:
    input:
        bam = "pipeline/alignment_pooled/pooled.{aligner}.bam"
    output:
        expand("pipeline/alignment_pooled/pooled.subsampled.{fraction}.{{aligner}}.bam", fraction=range(20, 100, 30))
    threads: 10
    resources:
        mem_mb = 400000,
        time_min = 1000,
        io_gb = 100
    params:
        tmpdir = "1000",
        outdir = "pipeline/alignment_pooled/"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bash workflow/scripts/subsample.sh {input.bam} 20 90 30 {threads} {params.outdir}"
