#!/bin/bash

function usage {
  echo "Usage: $0 [--indir] [--outdir] [--suffix] [--threads]"
  echo
  echo "Runs STAR alignment program on trimmed fastq files in a specified input directory and writes .sam files to an output directory. Uses hg38 reference genome by default"
  echo
  echo "  -i/--indir    input directory of trimmed fastq files"
  echo "  -o/--outdir   output directory for sam files"
  echo "  -s/--suffix   suffix denoting forward paired reads in the input directory, e.g. _1P.fastq"
  echo "  -t/--threads  number of threads"
  echo "  Version control: Nick Geraci July 27, 2017"
  echo
  exit 1
}

if [ "$1" == "-h" ]
then
  usage
fi

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
  -i|--indir) # directory for trimmed files
  INDIR="$2"
  shift
  ;;
  -o|--outdir) # output directory
  OUTDIR="$2"
  shift
  ;;
  -t|--threads) # number of threads
  THREADS="$2"
  shift
  ;;
  *)
  
  ;;
esac
shift
done

GENOMEDIR=/mnt/disks/disk-2/Tools/STAR-master/STARindex/hg38
cd $INDIR
INFILES=SRR*.fastq

# loop over files to align
for f in $INFILES
do
  OUTFILE=${f%%.fastq}"_STAR"
  OUTPATH="$OUTDIR/$OUTFILE"
  echo
  echo ALIGNING $f
  echo
  '/mnt/disks/disk-2/Tools/STAR-master/bin/Linux_x86_64/STAR' --genomeDir $GENOMEDIR --readFilesIn $f --outFileNamePrefix $OUTPATH --runThreadN $THREADS
done
