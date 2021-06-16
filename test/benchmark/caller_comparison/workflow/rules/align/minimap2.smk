localrules: sort_minimap2_alignments

#align reads (MD tags)
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
        "../../../../envs/minimap2.yaml"
    shell:
        "minimap2 -ax {params.preset} {params.options} -t {threads} --MD -Y {input.genome} {input.fq} > {output}"

#sort alignments
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
        "../../../../envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"
