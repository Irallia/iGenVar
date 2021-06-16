#!/bin/bash
set -euo pipefail

BAM=$1
BAI=$2
mkdir $MXQ_JOB_TMPDIR/input
echo "Copy $BAM to $MXQ_JOB_TMPDIR/input"
cp $BAM $MXQ_JOB_TMPDIR/input
cp $BAI $MXQ_JOB_TMPDIR/input

FROM=$3
TO=$4
STEP=$5
MINLEN=$6
THREADS=$7
BAM_NAME="$(basename -- $BAM)"
mkdir -p $MXQ_JOB_TMPDIR/output
echo "Start jobs"
seq $FROM $STEP $TO | xargs -I MINSUPP -P $THREADS -n 1 sniffles --mapped_reads $MXQ_JOB_TMPDIR/input/$BAM_NAME --min_length $MINLEN --min_support MINSUPP --vcf $MXQ_JOB_TMPDIR/output/raw_MINSUPP.vcf --threads 1 --genotype

OUTDIR=$8
echo "Copy back to $OUTDIR"
seq $FROM $STEP $TO | xargs -I MINSUPP -n 1 cp $MXQ_JOB_TMPDIR/output/raw_MINSUPP.vcf $OUTDIR
