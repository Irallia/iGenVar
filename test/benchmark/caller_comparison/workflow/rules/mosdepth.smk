
rule mosdepth_get:
    input:
        bam = "pipeline/alignment_pooled/{subsample}.{aligner}.bam",
        bai = "pipeline/alignment_pooled/{subsample}.{aligner}.bam.bai"
    threads: 4
    output:
        protected("pipeline/mosdepth/{subsample}.{aligner}.mosdepth.global.dist.txt"),
        protected("pipeline/mosdepth/{subsample}.{aligner}.regions.bed.gz"),
    resources:
        io_gb = 100
    conda:
        "../envs/mosdepth.yaml"
    params:
        windowsize = 500,
        workdir = "pipeline/mosdepth/{subsample}.{aligner}"
    log:
        "pipeline/logs/mosdepth/mosdepth_{subsample}.{aligner}.log"
    shell:
        "mosdepth --threads {threads} \
                  -n \
                  --by {params.windowsize} \
                  {params.workdir} {input.bam} 2> {log}"


# rule mosdepth_combine:
#     input:
#         expand("{{aligner}}/mosdepth/{subsample}.regions.bed.gz", subsample=["pooled.subsampled.{0}".format(cov) for cov in range(10, 91, 10)] + ["pooled"])
#     output:
#         "{aligner}/mosdepth/regions.combined.gz"
#     log:
#         "logs/{aligner}/mosdepth/mosdepth_combine.log"
#     shell:
#         os.path.join(workflow.basedir, "scripts/combine_mosdepth.py") + \
#             " {input} -o {output} 2> {log}"


# rule mosdepth_global_plot:
#     input:
#         expand("{{aligner}}/mosdepth/{subsample}.mosdepth.global.dist.txt", subsample=["pooled.subsampled.{0}".format(cov) for cov in range(10, 91, 10)] + ["pooled"])
#     output:
#         "{aligner}/mosdepth_global_plot/global.html"
#     log:
#         "logs/{aligner}/mosdepth/mosdepth_global_plot.log"
#     shell:
#         os.path.join(workflow.basedir, "scripts/mosdepth_plot-dist.py") + \
#             " {input} -o {output} 2> {log}"


rule compute_mean_converage:
    input:
        expand("pipeline/mosdepth/{subsample}.{{aligner}}.mosdepth.global.dist.txt", subsample=SUBSAMPLES)
    output:
        "pipeline/mosdepth/mean_coverages.{aligner}.txt"
    log:
        "pipeline/logs/mosdepth/mosdepth_mean_coverage.{aligner}.log"
    run:
        shell("rm -f {output}")
        for f in input:
            short_name = f.split("/")[2][:-25]
            shell("echo -n -e \"{short_name}\t\" >> {output}")
            shell("python workflow/scripts/mosdepth_mean-depth.py {f} >> {output}")