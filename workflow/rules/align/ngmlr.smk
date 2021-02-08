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
        "../../envs/ngmlr.yaml"
    shell:
        "cat {input.fq} | \
         ngmlr -presets {params.preset} -t {threads} -r {input.genome} | \
         samtools sort -@ {threads} -o {output} -"