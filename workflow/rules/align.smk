localrules: pbmm_index, index_alignment, alignment_stats, pool_samples

def get_samples(wildcards):
    return config["samples"][wildcards.sample]


rule run_alignments_minimap2:
    input:
        fq = get_samples,
        genome = config["reference"]
    output:
        temp("pipeline/alignments/{sample}.minimap2.unsorted.bam")
    resources:
        mem_mb = 50000,
        time_min = 1500,
        io_gb = 100
    params:
        preset = config["parameters"]["minimap_preset"],
        options = config["parameters"]["minimap_options"]
    threads: 10
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax {params.preset} {params.options} -t {threads} --MD -Y {input.genome} {input.fq} > {output}"

rule sort_minimap2_alignments:
    input:
        "pipeline/alignments/{sample}.minimap2.unsorted.bam"
    output:
        "pipeline/alignments/{sample}.minimap2.bam"
    resources:
        mem_mb = 50000,
        time_min = 1500,
        io_gb = 100
    threads: 10
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule pbmm_index:
    input:
        genome = config["reference"]
    output:
        index = config["reference"] + ".mmi"
    params:
        preset = config["parameters"]["pbmm_preset"]
    resources:
        io_gb = 100
    threads: 2
    conda:
        "../envs/pbmm2.yaml"
    shell:
        "pbmm2 index --num-threads {threads} --preset {params.preset} \
        {input.genome} {output.index}"

rule run_alignments_pbmm2:
    input:
        fq = get_samples,
        index = config["reference"] + ".mmi"
    output:
        bam = temp("pipeline/alignments/{sample}.pbmm2.bam")
    threads: 10
    resources:
        mem_mb = 50000,
        time_min = 1500,
        io_gb = 100
    params:
        sample = "{sample}",
        preset = config["parameters"]["pbmm_preset"]
    conda:
        "../envs/pbmm2.yaml"
    shell:
        """
        pbmm2 align --preset {params.preset} -j {threads} \
        --sort --rg '@RG\tID:rg1a\tSM:{params.sample}' --sample HG2 \
        {input.index} {input.fq} {output.bam}
        """

rule run_alignments_ngmlr:
    input:
        fq = get_samples,
        genome = config["reference"]
    output:
        temp("pipeline/alignments/{sample}.ngmlr.bam")
    resources:
        mem_mb = 50000,
        time_min = 1500,
        io_gb = 100
    params:
        preset = config["parameters"]["ngmlr_preset"]
    threads: 10
    conda:
        "../envs/ngmlr.yaml"
    shell:
        "zcat {input.fq} | \
         ngmlr --presets {params.preset} -t {threads} -r {input.genome} | \
         samtools sort -@ {threads} -o {output} -"

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
        "samtools merge -r {output} {input}"

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
        tmpdir = "700",
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
        tmpdir = "700",
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
        tmpdir = "700",
        outdir = "pipeline/alignment_pooled/"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bash workflow/scripts/subsample.sh {input.bam} 20 90 30 {threads} {params.outdir}"
