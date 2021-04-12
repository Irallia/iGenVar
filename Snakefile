import os
import gzip
import sys
import shutil

configfile: "config/config.yaml"

ALIGNERS=["minimap2", "pbmm2"]
SUBSAMPLES=["pooled"] + [("pooled.subsampled." + str(fraction)) for fraction in range(10, 100, 10)]
VCFS=["giab", "giab.gt", "giab.seq", "giab.gt.seq"]
SVIM_THRESHOLDS = [1, 3, 4, 5, 6, 7, 8, 9, 10, 14, 18, 22, 26, 30, 40, 50, 60]

wildcard_constraints:
    aligner="minimap2|ngmlr|pbmm2|lra"

include: "workflow/rules/align.smk"
include: "workflow/rules/mosdepth.smk"
include: "workflow/rules/callers.smk"
include: "workflow/rules/eval.smk"
include: "workflow/rules/plots.smk"

##### Target rules #####

rule all:
    input:
        #Alignments
        expand("pipeline/alignment_stats/alignment_stats.{aligner}.txt", aligner=ALIGNERS),
        expand("pipeline/mosdepth/mean_coverages.{aligner}.txt", aligner=ALIGNERS),
        #SV lengths
        expand("pipeline/SV-plots/{aligner}/pooled/SV-length_pbsv_{minscore}.png", aligner=ALIGNERS, minscore=[3, 5, 7]),
        expand("pipeline/SV-plots/{aligner}/pooled/SV-length_Sniffles_{minscore}.png", aligner=ALIGNERS, minscore=[3, 5, 7]),
        expand("pipeline/SV-plots/{aligner}/pooled/SV-length_SVIM_1000_900_1.0_0.5_{minscore}.png", aligner=ALIGNERS, minscore=[3, 5, 7]),
        expand("pipeline/SV-plots/{aligner}/pooled/SV-counts.merged.png", aligner=ALIGNERS),
        #Evaluation
        expand("pipeline/eval/{aligner}/results.{aligner}.all.png", aligner=ALIGNERS),
        expand("pipeline/eval/{aligner}/results.{aligner}.tools.{vcf}.png", aligner=ALIGNERS, vcf=VCFS),
        expand("pipeline/eval/{aligner}/results.{aligner}.coverages.{vcf}.png", aligner=ALIGNERS, vcf=VCFS),
        expand("pipeline/eval/{aligner}/results.{aligner}.svim.parameters.png", aligner=["pbmm2"]),
