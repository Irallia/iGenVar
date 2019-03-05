rule sort_calls:
    input:
        "{aligner}/{caller}_calls/{sample}.min_{minscore}.vcf"
    output:
        "{aligner}/{caller,sniffles|pbsv}_calls/{sample}.min_{minscore,[0-9]+}.sorted.vcf"
    threads: 1
    log:
        "logs/{aligner}/bcftools_sort/sorting_{caller}_{sample}_{minscore}.log"
    shell:
        "bcftools sort {input} > {output} 2> {log}"

rule sort_calls_svim:
    input:
        "{aligner}/svim_calls/{sample}/{parameters}/min_{minscore}.truvari.vcf"
    output:
        "{aligner}/svim_calls/{sample}/{parameters}/min_{minscore,[0-9]+}.truvari.sorted.vcf"
    threads: 1
    log:
        "logs/{aligner}/bcftools_sort/sorting_svim_{sample}_{parameters}_{minscore}.truvari.log"
    shell:
        "bcftools sort {input} > {output} 2> {log}"

rule bgzip:
    input:
        "{name}.sorted.vcf"
    output:
        "{name}.sorted.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"


rule tabix:
    input:
        "{name}.sorted.vcf.gz"
    output:
        "{name}.sorted.vcf.gz.tbi"
    shell:
        "tabix {input}"


rule callset_eval_svim:
    input:
        genome = config["genome"],
        truth_vcf = config["truth"]["vcf"],
        truth_bed = config["truth"]["bed"],
        calls = "{aligner}/svim_calls/{sample}/{parameters}/min_{minscore}.truvari.sorted.vcf.gz",
        index = "{aligner}/svim_calls/{sample}/{parameters}/min_{minscore}.truvari.sorted.vcf.gz.tbi"
    output:
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/summary.txt"
    params:
        out_dir="{aligner}/svim_results/{sample}/{parameters}/{minscore}"
    threads: 1
    log:
        "logs/{aligner}/truvari/pooled.svim.{sample}.{parameters}.{minscore}.log"
    shell:
        "rm -r {params.out_dir} && /project/pacbiosv/bin/truvari/truvari.py -f {input.genome}\
                    -b {input.truth_vcf} -c {input.calls} -o {params.out_dir}\
                    --passonly --includebed {input.truth_bed} --giabreport -r 1000 -p 0.00 2> {log}"

rule callset_eval_sniffles:
    input:
        genome = config["genome"],
        truth_vcf = config["truth"]["vcf"],
        truth_bed = config["truth"]["bed"],
        calls = "{aligner}/sniffles_calls/{sample}.min_{minscore}.sorted.vcf.gz",
        index = "{aligner}/sniffles_calls/{sample}.min_{minscore}.sorted.vcf.gz.tbi"
    output:
        "{aligner}/sniffles_results/{sample}/{minscore}/summary.txt"
    params:
        out_dir="{aligner}/sniffles_results/{sample}/{minscore}"
    threads: 1
    log:
        "logs/{aligner}/truvari/pooled.sniffles.{sample}.{minscore}.log"
    shell:
        "rm -r {params.out_dir} && /project/pacbiosv/bin/truvari/truvari.py -f {input.genome}\
                    -b {input.truth_vcf} -c {input.calls} -o {params.out_dir}\
                    --passonly --includebed {input.truth_bed} --giabreport -r 1000 -p 0.00 2> {log}"

rule callset_eval_pbsv:
    input:
        genome = config["genome"],
        truth_vcf = config["truth"]["vcf"],
        truth_bed = config["truth"]["bed"],
        calls = "{aligner}/pbsv_calls/{sample}.min_{minscore}.sorted.vcf.gz",
        index = "{aligner}/pbsv_calls/{sample}.min_{minscore}.sorted.vcf.gz.tbi"
    output:
        "{aligner}/pbsv_results/{sample}/{minscore}/summary.txt"
    params:
        out_dir="{aligner}/pbsv_results/{sample}/{minscore}"
    threads: 1
    log:
        "logs/{aligner}/truvari/pooled.pbsv.{sample}.{minscore}.log"
    shell:
        "rm -r {params.out_dir} && /project/pacbiosv/bin/truvari/truvari.py -f {input.genome}\
                    -b {input.truth_vcf} -c {input.calls} -o {params.out_dir}\
                    --passonly --includebed {input.truth_bed} --giabreport -r 1000 -p 0.00 2> {log}"

rule reformat_truvari_results:
    input:
        "{aligner}/{caller}_results/{sample}/{minscore}/summary.txt"
    output:
        "{aligner}/{caller,sniffles|pbsv}_results/{sample}/{minscore}/pr_rec.txt"
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", \"{wildcards.sample}\", {wildcards.minscore}, $1, $2 }}' > {output}"

rule reformat_truvari_results_svim:
    input:
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/summary.txt"
    output:
        "{aligner}/svim_results/{sample}/{parameters}/{minscore}/pr_rec.txt"
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"svim\", \"{wildcards.sample}\", {wildcards.minscore}, $1, $2 }}' > {output}"


rule cat_truvari_results_all:
    input:
        expand("{{aligner}}/svim_results/{{sample}}/{minscore}/pr_rec.txt", minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"])),
        expand("{{aligner}}/sniffles_results/{{sample}}/{minscore}/pr_rec.txt", minscore=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"])),
        expand("{{aligner}}/pbsv_results/{{sample}}/{minscore}/pr_rec.txt", minscore=range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"], config["minimums"]["pbsv_step"]))
    output:
        "{aligner}/eval/{sample}/all_results.txt"
    threads: 1
    shell:
        "cat {input} > {output}"

rule cat_truvari_results_nosniffles:
    input:
        expand("{{aligner}}/svim_results/{{sample}}/{minscore}/pr_rec.txt", minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"])),
        expand("{{aligner}}/pbsv_results/{{sample}}/{minscore}/pr_rec.txt", minscore=range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"], config["minimums"]["pbsv_step"]))
    output:
        "{aligner}/eval/{sample}/svim_pbsv_results.txt"
    threads: 1
    shell:
        "cat {input} > {output}"

rule cat_truvari_results_nopbsv:
    input:
        expand("{{aligner}}/svim_results/{{sample}}/{{parameters}}/{minscore}/pr_rec.txt", minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"])),
        expand("{{aligner}}/sniffles_results/{{sample}}/{minscore}/pr_rec.txt", minscore=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"]))
    output:
        "{aligner}/eval/{sample}/{parameters}/svim_sniffles_results.txt"
    threads: 1
    shell:
        "cat {input} > {output}"

rule cat_truvari_results_svim_multiple_coverages:
    input:
        expand("{{aligner}}/svim_results/{{sample}}.subsampled.{fraction}/{{parameters}}/{minscore}/pr_rec.txt", fraction=range(10, 91, 10), minscore=range(config["minimums"]["svim_from"], config["minimums"]["svim_to"], config["minimums"]["svim_step"]))
    output:
        "{aligner}/eval/{sample}/{parameters}/svim_results_multiple_coverages.txt"
    threads: 1
    shell:
        "cat {input} > {output}"

rule cat_truvari_results_sniffles_multiple_coverages:
    input:
        expand("{{aligner}}/sniffles_results/{{sample}}.subsampled.{fraction}/{minscore}/pr_rec.txt", fraction=range(10, 91, 10), minscore=range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"], config["minimums"]["sniffles_step"]))
    output:
        "{aligner}/eval/{sample}/sniffles_results_multiple_coverages.txt"
    threads: 1
    shell:
        "cat {input} > {output}"

rule cat_truvari_results_pbsv_multiple_coverages:
    input:
        expand("{{aligner}}/pbsv_results/{{sample}}.subsampled.{fraction}/{minscore}/pr_rec.txt", fraction=range(10, 91, 10), minscore=range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"], config["minimums"]["pbsv_step"]))
    output:
        "{aligner}/eval/{sample}/pbsv_results_multiple_coverages.txt"
    threads: 1
    shell:
        "cat {input} > {output}"

rule plot_pr_tools:
    input:
        "{aligner}/eval/{sample}/all_results.txt"
    output:
        "{aligner}/eval/{sample}/tools_pr_all.png"
    threads: 1
    log:
        "logs/{aligner}/rplot/{sample}.tools.pr.log"
    shell:
        "Rscript --vanilla scripts/plot-pr-tools.R {input} {output} > {log}"

rule plot_pr_tools_nosniffles:
    input:
        "{aligner}/eval/{sample}/svim_pbsv_results.txt"
    output:
        "{aligner}/eval/{sample}/tools_pr_svim_pbsv.png"
    threads: 1
    log:
        "logs/{aligner}/rplot/{sample}.tools.pr.log"
    shell:
        "Rscript --vanilla scripts/plot-pr-tools.R {input} {output} > {log}"

rule plot_pr_tools_nopbsv:
    input:
        "{aligner}/eval/{sample}/{parameters}/svim_sniffles_results.txt"
    output:
        "{aligner}/eval/{sample}/{parameters}/tools_pr_svim_sniffles.png"
    threads: 1
    log:
        "logs/{aligner}/rplot/{sample}.{parameters}.tools.pr.log"
    shell:
        "Rscript --vanilla scripts/plot-pr-tools.R {input} {output} > {log}"

rule plot_pr_tools_multiple_coverages:
    input:
        "{aligner}/eval/{sample}/{caller}_results_multiple_coverages.txt"
    output:
        "{aligner}/eval/{sample}/{caller,sniffles|pbsv}_pr_multiple_coverages.png"
    threads: 1
    log:
        "logs/{aligner}/rplot/{sample}.{caller}.pr.coverages.log"
    shell:
        "Rscript --vanilla scripts/plot-pr-tools-coverages.R {input} {output} > {log}"

rule plot_pr_tools_multiple_coverages_svim:
    input:
        "{aligner}/eval/{sample}/{parameters}/svim_results_multiple_coverages.txt"
    output:
        "{aligner}/eval/{sample}/{parameters}/svim_pr_multiple_coverages.png"
    threads: 1
    log:
        "logs/{aligner}/rplot/{sample}.{parameters}.svim.pr.coverages.log"
    shell:
        "Rscript --vanilla scripts/plot-pr-tools-coverages.R {input} {output} > {log}"
