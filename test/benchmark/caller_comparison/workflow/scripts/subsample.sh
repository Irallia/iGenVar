#!/bin/bash
set -euo pipefail

BAM=$1
mkdir -p $MXQ_JOB_TMPDIR/input
echo "Copy $BAM to $MXQ_JOB_TMPDIR/input"
cp $BAM $MXQ_JOB_TMPDIR/input

FROM=$2
TO=$3
STEP=$4
THREADS=$5
OUTDIR=$6
BAM_NAME="$(basename -- $BAM)"
A="$(cut -d'.' -f1 <<< "$BAM_NAME")"
B="$(cut -d'.' -f2 <<< "$BAM_NAME")"
C="$(cut -d'.' -f3 <<< "$BAM_NAME")"
mkdir -p $MXQ_JOB_TMPDIR/output
echo "Start jobs"
seq $FROM $STEP $TO | xargs -I FRAC -P $THREADS -n 1 samtools view -s 10.FRAC -b $MXQ_JOB_TMPDIR/input/$BAM_NAME -o $MXQ_JOB_TMPDIR/output/$A.subsampled.FRAC.$B.$C

echo "Move back to $OUTDIR"
seq $FROM $STEP $TO | xargs -I FRAC -n 1 -P 1 mv $MXQ_JOB_TMPDIR/output/$A.subsampled.FRAC.$B.$C $OUTDIR
