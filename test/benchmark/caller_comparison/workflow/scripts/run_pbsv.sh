#!/bin/bash
set -euo pipefail

BAM=$1
BAI=$2
GENOME=$3
mkdir $MXQ_JOB_TMPDIR/input
echo "Copy $BAM and $GENOME to $MXQ_JOB_TMPDIR/input"
cp $BAM $MXQ_JOB_TMPDIR/input
cp $BAI $MXQ_JOB_TMPDIR/input
cp $GENOME $MXQ_JOB_TMPDIR/input

FROM=$4
TO=$5
STEP=$6
MINLEN=$7
THREADS=$8
BAM_NAME="$(basename -- $BAM)"
GENOME_NAME="$(basename -- $GENOME)"
mkdir -p $MXQ_JOB_TMPDIR/output
echo "Start pbsv discover"
pbsv discover $MXQ_JOB_TMPDIR/input/$BAM_NAME $MXQ_JOB_TMPDIR/output/signatures.svsig.gz
echo "Start pbsv call"
seq $FROM $STEP $TO | xargs -I MINSUPP -P $THREADS -n 1 pbsv call -t INS,DEL -j 1 --min-sv-length $MINLEN --max-ins-length 100K --call-min-reads-one-sample MINSUPP --call-min-reads-all-samples MINSUPP --call-min-reads-per-strand-all-samples 0 --call-min-bnd-reads-all-samples 0 --call-min-read-perc-one-sample 0 $MXQ_JOB_TMPDIR/input/$GENOME_NAME $MXQ_JOB_TMPDIR/output/signatures.svsig.gz $MXQ_JOB_TMPDIR/output/min_MINSUPP.vcf

OUTDIR=$9
echo "Copy back to $OUTDIR"
seq $FROM $STEP $TO | xargs -I MINSUPP -n 1 cp $MXQ_JOB_TMPDIR/output/min_MINSUPP.vcf $OUTDIR

