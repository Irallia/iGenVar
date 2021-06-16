localrules: plot_pr_all_results, plot_pr_tools, plot_pr_coverages, plot_pr_svim_parameters, run_pbsv1_all_types, run_pbsv2_all_types, SV_length_plot_pbsv, SV_length_plot_sniffles, SV_length_plot_svim, merge_counts, plot_counts

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
        "pipeline/eval/{aligner}/full_results.txt"
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
        png = "pipeline/eval/{aligner}/results.{aligner}.coverages.{vcf}.png",
        tsv = "pipeline/eval/{aligner}/results.{aligner}.coverages.{vcf}.tsv"
    threads: 1
    log:
        "pipeline/logs/rplot.coverages.{aligner}.{vcf}.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_coverages.R {input} {wildcards.vcf} {output.png} {output.tsv} > {log}"

rule plot_pr_svim_parameters:
    input:
        "pipeline/eval/{aligner}/svim_parameter_results.txt"
    output:
        "pipeline/eval/{aligner}/results.{aligner}.svim.parameters.png"
    threads: 1
    log:
        "pipeline/logs/rplot.all.{aligner}.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_svim_parameters.R {input} {output} > {log}"

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

#run pbsv for all SV types
rule run_pbsv1_all_types:
    input:
        bam = "pipeline/alignment_pooled/{data}.pbmm2.bam",
        bai = "pipeline/alignment_pooled/{data}.pbmm2.bam.bai",
        genome = config["reference"],
    output:
        svsig = temp("pipeline/pbsv/{aligner}/{data}/svsig_all_types.svsig.gz")
    resources:
        mem_mb = 10000,
        time_min = 600,
        io_gb = 100
    threads: 1
    conda:
        "../envs/pbsv.yaml"
    shell:
        """
        pbsv discover {input.bam} {output.svsig}
        """

rule run_pbsv2_all_types:
    input:
        svsig = "pipeline/pbsv/{aligner}/{data}/svsig_all_types.svsig.gz",
        genome = config["reference"],
    output:
        vcf = "pipeline/pbsv/{aligner}/{data}/min_{minsupport,[0-9]+}_all_types.vcf"
    params:
        min_sv_size = config["parameters"]["min_sv_size"]
    resources:
        mem_mb = 10000,
        time_min = 600,
        io_gb = 100
    threads: 1
    conda:
        "../envs/pbsv.yaml"
    shell:
        """
        pbsv call -j 1 \
        --min-sv-length {params.min_sv_size} \
        --max-ins-length 100K \
        --call-min-reads-one-sample {wildcards.minsupport} \
        --call-min-reads-all-samples {wildcards.minsupport} \
        --call-min-reads-per-strand-all-samples 0 \
        --call-min-bnd-reads-all-samples 0 \
        --call-min-read-perc-one-sample 0 {input.genome} {input.svsig} {output.vcf}
        """

rule SV_length_plot_pbsv:
    input:
        "pipeline/pbsv/{aligner}/{data}/min_{minimum}_all_types.vcf"
    output:
        plot = "pipeline/SV-plots/{aligner}/{data}/SV-length_pbsv_{minimum}.png",
        counts = "pipeline/SV-plots/{aligner}/{data}/SV-counts_pbsv_{minimum}.txt",
    log:
        "logs/svplot/svlength_pbsv_{aligner}_{data}_{minimum}.log"
    conda:
        "../envs/cyvcf2.yaml"
    shell:
        "python workflow/scripts/SV-length-plot.py {input} --output {output.plot} --counts {output.counts} --filter 'hs37d5' --tool pbsv 2> {log}"

rule SV_length_plot_sniffles:
    input:
        "pipeline/Sniffles/{aligner}/{data}/min_{minimum}.vcf"
    output:
        plot = "pipeline/SV-plots/{aligner}/{data}/SV-length_Sniffles_{minimum}.png",
        counts = "pipeline/SV-plots/{aligner}/{data}/SV-counts_Sniffles_{minimum}.txt",
    log:
        "logs/svplot/svlength_Sniffles_{aligner}_{data}_{minimum}.log"
    conda:
        "../envs/cyvcf2.yaml"
    shell:
        "python workflow/scripts/SV-length-plot.py {input} --output {output.plot} --counts {output.counts} --filter 'hs37d5'  --tool Sniffles 2> {log}"

#run SVIM without converting DUP to INS
rule run_svim_all_types:
    input:
        genome = config["reference"],
        bam="pipeline/alignment_pooled/{data}.{aligner}.bam",
        bai="pipeline/alignment_pooled/{data}.{aligner}.bam.bai"
    output:
        "pipeline/SVIM/{aligner}/{data}/{pmd,[0-9]+}_{pdn,[0-9]+}_{edn,[0-9\.]+}_{cmd,[0-9\.]+}_all_types/variants.vcf"
    resources:
        mem_mb = 10000,
        time_min = 600,
        io_gb = 100
    params:
        working_dir = "pipeline/SVIM/{aligner}/{data}/{pmd}_{pdn}_{edn}_{cmd}_all_types/",
        min_sv_size = config["parameters"]["min_sv_size"]
    threads: 1
    shell:
        "/home/heller_d/bin/anaconda3/envs/svim_test4/bin/python /home/heller_d/bin/anaconda3/envs/svim_test4/bin/svim alignment --sample {wildcards.data} \
         --partition_max_distance {wildcards.pmd} \
         --position_distance_normalizer {wildcards.pdn} \
         --edit_distance_normalizer {wildcards.edn} \
         --cluster_max_distance {wildcards.cmd} \
         --min_sv_size {params.min_sv_size} \
         --segment_gap_tolerance 20 \
         --segment_overlap_tolerance 20 \
         --read_names \
         --max_sv_size 1000000 \
         --verbose \
         {params.working_dir} {input.bam} {input.genome}"

rule SV_length_plot_svim:
    input:
        "pipeline/SVIM/{aligner}/{data}/{parameters}_all_types/variants.vcf"
    output:
        plot = "pipeline/SV-plots/{aligner}/{data}/SV-length_SVIM_{parameters}_{minimum}.png",
        counts = "pipeline/SV-plots/{aligner}/{data}/SV-counts_SVIM_{parameters}_{minimum}.txt",
    log:
        "logs/svplot/svlength_SVIM_{aligner}_{data}_{parameters}_{minimum}.log"
    conda:
        "../envs/cyvcf2.yaml"
    shell:
        "python workflow/scripts/SV-length-plot.py {input} --min_score {wildcards.minimum} --output {output.plot} --counts {output.counts} --filter 'hs37d5'  --tool SVIM 2> {log}"

rule merge_counts:
    input:
        svim = "pipeline/SV-plots/{aligner}/{data}/SV-counts_SVIM_1000_900_1.0_0.5_7.txt",
        sniffles = "pipeline/SV-plots/{aligner}/{data}/SV-counts_Sniffles_5.txt",
        pbsv = "pipeline/SV-plots/{aligner}/{data}/SV-counts_pbsv_5.txt",
    output:
        "pipeline/SV-plots/{aligner}/{data}/SV-counts.merged.txt"
    shell:
        "cat {input} | grep -v '#' > {output}"

rule plot_counts:
    input:
        "pipeline/SV-plots/{aligner}/{data}/SV-counts.merged.txt"
    output:
        png = "pipeline/SV-plots/{aligner}/{data}/SV-counts.merged.png",
        tsv = "pipeline/SV-plots/{aligner}/{data}/SV-counts.merged.tsv"
    shell:
        "Rscript --vanilla workflow/scripts/plot_counts.R {input} {output.png} {output.tsv}"
