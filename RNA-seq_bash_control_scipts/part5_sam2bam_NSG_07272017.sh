#!/bin/bash

function usage {
  echo "Usage: $0 [--indir] [--outdir] [--suffix] [--threads] [--memperthread]"
  echo
  echo "Converts .sam files in a specified input directory to sorted .bam files in an output directory."
  echo
  echo "  -i/--indir          input directory of sam files"
  echo "  -o/--outdir         output directory for sorted bam files"
  echo "  -s/--suffix         suffix denoting STAR output sam files for conversion, e.g. _STARAligned.out.sam"
  echo "  -t/--threads        number of threads"
  echo "  -m/--memperthread   memory per thread, e.g. 2G"
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
  -i|--indir) # input directory
  INDIR="$2"
  shift
  ;;
  -o|--outdir) # output directory
  OUTDIR="$2"
  shift
  ;;
  -s|--suffix) # _STARAligned.out.sam
  SUFFIX="$2"
  shift
  ;;
  -t|--threads) # number of threads
  THREADS="$2"
  shift
  ;;
  -m|--memperthread) # usually 2G for 2 gigabytes
  MEM="$2"
  shift
  ;;
  *)
  
  ;;
esac
shift
done

cd $INDIR
INFILES=*$SUFFIX

for f in $INFILES
do
  OUTFILE=${f%%$SUFFIX}"_sorted"
  OUTPATH="$OUTDIR/$OUTFILE"
  echo
  echo SORTING $f
  echo
  '/mnt/disks/disk-2/Tools/samtools-1.5/samtools' view -@ $THREADS -b -u -S $f | '/mnt/disks/disk-2/Tools/samtools-1.5/samtools' sort -o $OUTPATH -@ $THREADS -m $MEM
done
