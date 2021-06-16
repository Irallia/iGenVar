localrules: pbmm2_index, add_md_to_pbmm2_alignment

#index reference genome
rule pbmm2_index:
    input:
        genome = config["reference"]
    output:
        index = config["reference"] + ".pbmm2.mmi"
    params:
        preset = config["parameters"]["pbmm_preset"]
    resources:
        io_gb = 100
    threads: 2
    conda:
        "../../../../envs/pbmm2.yaml"
    shell:
        "pbmm2 index --num-threads {threads} --preset {params.preset} \
        {input.genome} {output.index}"

#align reads (sorted, read group)
rule run_alignments_pbmm2:
    input:
        fq = get_samples,
        index = ancient(config["reference"] + ".pbmm2.mmi")
    output:
        bam = temp("pipeline/alignments/{sample}.pbmm2.raw.bam")
    threads: 10
    resources:
        mem_mb = 50000,
        time_min = 1500,
        io_gb = 100
    params:
        sample = "{sample}",
        preset = config["parameters"]["pbmm_preset"]
    conda:
        "../../../../envs/pbmm2.yaml"
    shell:
        """
        pbmm2 align --preset {params.preset} -j {threads} \
        --sort --rg '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}' --sample HG2 \
        {input.index} {input.fq} {output.bam}
        """

#add MD tag
rule add_md_to_pbmm2_alignment:
    input:
        bam = "pipeline/alignments/{sample}.pbmm2.raw.bam",
        genome = config["reference"]
    output:
        bam = "pipeline/alignments/{sample}.pbmm2.bam"
    resources:
        io_gb = 100
    log:
        "pipeline/logs/calmd/calmd.pbmm2.{sample}.log"
    threads: 1
    shell:
        "samtools calmd -b {input.bam} {input.genome} 2> {log} > {output.bam}"
