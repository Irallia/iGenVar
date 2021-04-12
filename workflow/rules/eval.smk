localrules: bgzip, tabix, callset_eval_svim, callset_eval, reformat_truvari_results, reformat_truvari_results_svim, cat_truvari_results_all, cat_truvari_results_full, cat_truvari_results_svim_parameters

def get_vcf(wildcards):
    return config["truth"][wildcards.vcf.split(".")[0]]

rule bgzip:
    input:
        "{name}.vcf"
    output:
        "{name}.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"


rule tabix:
    input:
        "{name}.vcf.gz"
    output:
        "{name}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"


rule callset_eval_svim:
    input:
        genome = config["reference"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "pipeline/SVIM/{aligner}/{data}/{parameters}/min_{minscore}.indel.vcf.gz",
        index = "pipeline/SVIM/{aligner}/{data}/{parameters}/min_{minscore}.indel.vcf.gz.tbi"
    output:
        summary="pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}/summary.txt"
    params:
        out_dir="pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}",
        vcf="{vcf}"
    threads: 1
    log:
        log="logs/truvari/truvari.svim.{data}.{aligner}.{parameters}.{minscore}.{vcf}.log"
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/run_truvari.py"


rule callset_eval:
    input:
        genome = config["reference"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "pipeline/{caller}/{aligner}/{data}/min_{minscore}.indel.vcf.gz",
        index = "pipeline/{caller}/{aligner}/{data}/min_{minscore}.indel.vcf.gz.tbi"
    output:
        summary="pipeline/{caller,Sniffles|pbsv}_results/{aligner}/{data}/{minscore}/{vcf}/summary.txt"
    params:
        out_dir="pipeline/{caller}_results/{aligner}/{data}/{minscore}/{vcf}",
        vcf="{vcf}"
    threads: 1
    log:
        log="logs/truvari/truvari.{caller}.{data}.{aligner}.{minscore}.{vcf}.log"
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/run_truvari.py"


rule reformat_truvari_results:
    input:
        "pipeline/{caller}_results/{aligner}/{data}/{minscore}/{vcf}/summary.txt"
    output:
        "pipeline/{caller,Sniffles|pbsv}_results/{aligner}/{data}/{minscore}/{vcf}/pr_rec.txt"
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", \"{wildcards.aligner}\", \"{wildcards.data}\", \"{wildcards.vcf}\", {wildcards.minscore}, $1, $2 }}' > {output}"

rule reformat_truvari_results_svim:
    input:
        "pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}/summary.txt"
    output:
        "pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}/pr_rec.txt"
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"SVIM\", \"{wildcards.aligner}\", \"{wildcards.data}\", \"{wildcards.vcf}\", {wildcards.minscore}, $1, $2 }}' > {output}"


rule cat_truvari_results_all:
    input:
        svim = expand("pipeline/SVIM_results/{{aligner}}/{data}/1000_900_0.3/{minscore}/{vcf}/pr_rec.txt", 
                          data = SUBSAMPLES, 
                          minscore=[0] + list(range(1, 60, 2)), vcf=VCFS),
        sniffles = expand("pipeline/Sniffles_results/{{aligner}}/{data}/{minscore}/{vcf}/pr_rec.txt", 
                          data = SUBSAMPLES, 
                          minscore=list(range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"]+1, config["minimums"]["sniffles_step"])),
                          vcf=VCFS),
        pbsv = expand("pipeline/pbsv_results/{{aligner}}/{data}/{minscore}/{vcf}/pr_rec.txt", 
                          data = SUBSAMPLES, 
                          minscore=list(range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"]+1, config["minimums"]["pbsv_step"])), 
                          vcf=VCFS)
    output:
        svim = temp("pipeline/eval/{aligner}/svim.all_results.txt"),
        sniffles = temp("pipeline/eval/{aligner}/sniffles.all_results.txt"),
        pbsv = temp("pipeline/eval/{aligner}/pbsv.all_results.txt"),
        all = "pipeline/eval/{aligner}/all_results.txt"
    threads: 1
    run:
        shell("cat {input.svim} > {output.svim}")
        shell("cat {input.sniffles} > {output.sniffles}")
        shell("cat {input.pbsv} > {output.pbsv}")
        shell("cat {output.svim} {output.sniffles} {output.pbsv} > {output.all}")

rule cat_truvari_results_full:
    input:
        svim = expand("pipeline/SVIM_results/{{aligner}}/pooled/1000_900_0.3/{minscore}/{vcf}/pr_rec.txt", 
                          minscore=[0] + list(range(1, 60, 2)), vcf=VCFS),
        sniffles = expand("pipeline/Sniffles_results/{{aligner}}/pooled/{minscore}/{vcf}/pr_rec.txt", 
                          minscore=list(range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"]+1, config["minimums"]["sniffles_step"])),
                          vcf=VCFS),
        pbsv = expand("pipeline/pbsv_results/{{aligner}}/pooled/{minscore}/{vcf}/pr_rec.txt", 
                          minscore=list(range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"]+1, config["minimums"]["pbsv_step"])), 
                          vcf=VCFS)
    output:
        svim = temp("pipeline/eval/{aligner}/svim.full_results.txt"),
        sniffles = temp("pipeline/eval/{aligner}/sniffles.full_results.txt"),
        pbsv = temp("pipeline/eval/{aligner}/pbsv.full_results.txt"),
        all = "pipeline/eval/{aligner}/full_results.txt"
    threads: 1
    run:
        shell("cat {input.svim} > {output.svim}")
        shell("cat {input.sniffles} > {output.sniffles}")
        shell("cat {input.pbsv} > {output.pbsv}")
        shell("cat {output.svim} {output.sniffles} {output.pbsv} > {output.all}")

rule cat_truvari_results_svim_parameters:
    input:
        svim = expand("pipeline/SVIM_results/{{aligner}}/{data}/{pmd}_{dn}_{cmd}/{minscore}/{vcf}/pr_rec.txt", 
                          data = ["pooled", "pooled.subsampled.50"],
                          pmd = [500, 1000, 5000],
                          dn = [900],
                          cmd = [0.2, 0.3, 0.4],
                          minscore=[0] + list(range(1, 60, 2)),
                          vcf=VCFS)
    output:
        all = "pipeline/eval/{aligner}/svim_parameter_results.txt"
    threads: 1
    run:
        with open(output.all, 'w') as output_file:
            for f in input.svim:
                parameters = f.split("/")[4]
                with open(f, 'r') as input_file:
                    for line in input_file:
                        print("%s\t%s" % (parameters, line), file=output_file)
