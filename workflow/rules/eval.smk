localrules: bgzip, tabix, callset_eval_svim, callset_eval, reformat_truvari_results, reformat_truvari_results_svim, cat_truvari_results

def get_vcf(wildcards):
    return config["truth"][wildcards.vcf]

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
        out_dir="pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}"
    threads: 1
    log:
        log="logs/truvari/truvari.svim.{data}.{aligner}.{parameters}.{minscore}.{vcf}.log"
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/run_truvari.py"

rule callset_eval_svim_gtcomp:
    input:
        genome = config["reference"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "pipeline/SVIM/{aligner}/{data}/{parameters}/min_{minscore}.indel.vcf.gz",
        index = "pipeline/SVIM/{aligner}/{data}/{parameters}/min_{minscore}.indel.vcf.gz.tbi"
    output:
        summary="pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}.gt/summary.txt"
    params:
        out_dir="pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}.gt"
    threads: 1
    log:
        log="logs/truvari/truvari.svim.{data}.{aligner}.{parameters}.{minscore}.{vcf}.gt.log"
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/run_truvari_gtcomp.py"

rule callset_eval:
    input:
        genome = config["reference"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "pipeline/{caller}/{aligner}/{data}/min_{minscore}.indel.vcf.gz",
        index = "pipeline/{caller}/{aligner}/{data}/min_{minscore}.indel.vcf.gz.tbi"
    output:
        "pipeline/{caller,Sniffles|pbsv}_results/{aligner}/{data}/{minscore}/{vcf}/summary.txt"
    params:
        out_dir="pipeline/{caller}_results/{aligner}/{data}/{minscore}/{vcf}"
    threads: 1
    log:
        log="logs/truvari/truvari.{caller}.{data}.{aligner}.{minscore}.{vcf}.log"
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/run_truvari.py"

rule callset_eval_gtcomp:
    input:
        genome = config["reference"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "pipeline/{caller}/{aligner}/{data}/min_{minscore}.indel.vcf.gz",
        index = "pipeline/{caller}/{aligner}/{data}/min_{minscore}.indel.vcf.gz.tbi"
    output:
        "pipeline/{caller,Sniffles|pbsv}_results/{aligner}/{data}/{minscore}/{vcf}.gt/summary.txt"
    params:
        out_dir="pipeline/{caller}_results/{aligner}/{data}/{minscore}/{vcf}.gt"
    threads: 1
    log:
        log="logs/truvari/truvari.{caller}.{data}.{aligner}.{minscore}.{vcf}.gt.log"
    conda:
        "../envs/truvari.yaml"
    script:
        "../scripts/run_truvari_gtcomp.py"

rule reformat_truvari_results:
    input:
        "pipeline/{caller}_results/{aligner}/{data}/{minscore}/{vcf}/summary.txt"
    output:
        temp("pipeline/{caller,Sniffles|pbsv}_results/{aligner}/{data}/{minscore}/{vcf}/pr_rec.txt")
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", \"{wildcards.aligner}\", \"{wildcards.data}\", \"{wildcards.vcf}\", {wildcards.minscore}, $1, $2 }}' > {output}"

rule reformat_truvari_results_svim:
    input:
        "pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}/summary.txt"
    output:
        temp("pipeline/SVIM_results/{aligner}/{data}/{parameters}/{minscore}/{vcf}/pr_rec.txt")
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"SVIM\", \"{wildcards.aligner}\", \"{wildcards.data}\", \"{wildcards.vcf}\", {wildcards.minscore}, $1, $2 }}' > {output}"


rule cat_truvari_results:
    input:
        svim = expand("pipeline/SVIM_results/{{aligner}}/{data}/0.3/{minscore}/{vcf}/pr_rec.txt", 
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
        svim = temp("pipeline/eval/{aligner}/svim.results.txt"),
        sniffles = temp("pipeline/eval/{aligner}/sniffles.results.txt"),
        pbsv = temp("pipeline/eval/{aligner}/pbsv.results.txt"),
        all = "pipeline/eval/{aligner}/all_results.txt"
    threads: 1
    run:
        shell("cat {input.svim} > {output.svim}")
        shell("cat {input.sniffles} > {output.sniffles}")
        shell("cat {input.pbsv} > {output.pbsv}")
        shell("cat {output.svim} {output.sniffles} {output.pbsv} > {output.all}")

rule plot_pr_all_results:
    input:
        "pipeline/eval/{aligner}/all_results.txt"
    output:
        "pipeline/eval/{aligner}/results.{aligner}.all.png"
    threads: 1
    log:
        "pipeline/logs/rplot.all.{aligner}.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_all_results.R {input} {output} > {log}"

rule plot_pr_tools:
    input:
        "pipeline/eval/{aligner}/all_results.txt"
    output:
        "pipeline/eval/{aligner}/results.{aligner}.tools.{vcf}.png"
    threads: 1
    log:
        "pipeline/logs/rplot.tools.{aligner}.{vcf}.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_tools.R {input} {wildcards.vcf} {output} > {log}"

rule plot_pr_coverages:
    input:
        "pipeline/eval/{aligner}/all_results.txt"
    output:
        "pipeline/eval/{aligner}/results.{aligner}.coverages.{vcf}.png"
    threads: 1
    log:
        "pipeline/logs/rplot.coverages.{aligner}.{vcf}.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_coverages.R {input} {wildcards.vcf} {output} > {log}"

# rule plot_pr_coverages_bar:
#     input:
#         "pipeline/eval/{aligner}/all_results.txt"
#     output:
#         "pipeline/eval/{aligner}/results.{aligner}.coverages.bar.png"
#     threads: 1
#     log:
#         "pipeline/logs/rplot.coveragesbar.{aligner}.log"
#     shell:
#         "Rscript --vanilla workflow/scripts/plot_coverages_bar.R {input} {output} > {log}"
