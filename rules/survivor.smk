rule survivor_combine_callers:
    input:
        expand("{{aligner}}/{caller}_calls/pooled.vcf",
               caller=["sniffles", "svim", "nanosv"])
    output:
        vcf = temp("{aligner}/pooled_combined/genotypes.vcf"),
        fofn = temp("{aligner}/pooled_combined/samples.fofn")
    params:
        distance = config["parameters"]["survivor_distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        "logs/{aligner}/all/survivor.log"
    shell:
        "ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"

rule sort_vcf:
    input:
        "{aligner}/pooled_combined/genotypes.vcf"
    output:
        "{aligner}/pooled_combined/genotypes.sorted.vcf"
    log:
        "logs/{aligner}/bcftools_sort/sorting_combined.log"
    threads: 8
    shell:
        "bcftools sort {input} > {output} 2> {log}"
