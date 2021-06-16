localrules: lra_index, add_md_to_lra_alignment, add_rg_to_lra_alignment

#index reference genome
rule lra_index:
    input:
        genome = config["reference"]
    output:
        mmi_index = config["reference"] + ".mmi",
        gli_index = config["reference"] + ".gli"
    params:
        preset = config["parameters"]["lra_preset"]
    resources:
        io_gb = 100
    threads: 1
    conda:
        "../../../../envs/lra.yaml"
    shell:
        "lra index {params.preset} \
        {input.genome}"

#align reads (sorted)
rule run_alignments_lra:
    input:
        fa = get_samples,
        genome = config["reference"],
        mmi_index = ancient(config["reference"] + ".mmi"),
        gli_index = ancient(config["reference"] + ".gli")
    output:
        bam = temp("pipeline/alignments/{sample}.lra.raw.bam")
    threads: 12
    resources:
        mem_mb = 50000,
        time_min = 1500,
        io_gb = 100
    params:
        preset = config["parameters"]["lra_preset"]
    conda:
        "../../../../envs/lra.yaml"
    shell:
        "lra align {params.preset} -p s -t 10 {input.genome} {input.fa} | \
        samtools sort -o {output.bam} -"

#add MD tag
rule add_md_to_lra_alignment:
    input:
        bam = "pipeline/alignments/{sample}.lra.raw.bam",
        genome = config["reference"]
    output:
        bam = temp("pipeline/alignments/{sample}.lra.md.bam")
    resources:
        io_gb = 100
    log:
        "pipeline/logs/calmd/calmd.lra.{sample}.log"
    threads: 1
    shell:
        "samtools calmd -b {input.bam} {input.genome} 2> {log} > {output.bam}"

#add read group
rule add_rg_to_lra_alignment:
    input:
        bam = "pipeline/alignments/{sample}.lra.md.bam"
    output:
        bam = "pipeline/alignments/{sample}.lra.bam"
    resources:
        io_gb = 100
    threads: 1
    shell:
        "samtools addreplacerg -r '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}' -o {output.bam} {input.bam}"
