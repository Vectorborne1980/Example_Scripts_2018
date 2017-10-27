#!/bin/bash

function usage {
  echo "Usage: $0 [--indir] [--outdir] [--suffix] [--threads] [--memperthread]"
  echo
  echo "Converts .sam files in a specified input directory to sorted .bam files in an output directory.  Most (?) efficient usage seems to be breaking the input files into 4 directories and creating 4 output directories, then running conversion in 4 terminals at once."
  echo
  echo "  -i/--indir          input directory of sam files"
  echo "  -o/--outdir         output directory for sorted bam files (need entire path)"
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

PATH="/usr/local/bin/anaconda3/bin:$PATH"

cd $INDIR
INFILES=*$SUFFIX

for f in $INFILES
do
  USORTFILE=${f%%$SUFFIX}"_unsorted"
  SORTFILE=${f%%$SUFFIX}"_sorted"
  USORTPATH="$OUTDIR/$USORTFILE"
  SORTPATH="$OUTDIR/$SORTFILE"
  echo
  echo VIEWING $f
  echo
  sambamba view -S -f bam -l 0 -t $THREADS -o $USORTPATH $f
  echo
  echo SORTING $USORTFILE
  echo
  sambamba sort -m $MEM --tmpdir $OUTDIR -o $SORTPATH -N -l 0 -p -t $THREADS $USORTPATH
  rm $USORTPATH
done
