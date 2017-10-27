#!/bin/bash

function usage {
  echo "Usage: $0 [--indir] [--outdir] [--threads]"
  echo
  echo "Runs parallel-fastq-dump on all .sra files in a specified input directory and writes .fastq files to an output directory. Breaks each file into blocks that are processed in parallel, and then iterates through all files. For whatever reason (PATH issues?), you may have to su before running this script. Also need parallel-fastq-dump installed for your username"
  echo
  echo "  -i/--indir    input directory of sra files"
  echo "  -o/--outdir   output directory for fastq files"
  echo "  -t/--threads  number of threads"
  echo "  Version control: Nick Geraci July 31, 2017"
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
  -i|--indir) # directory to read sra files from
  INDIR="$2"
  shift
  ;;
  -o|--outdir) # directory for output fastq files
  OUTDIR="$2"
  shift
  ;;
  -t|--threads) #number of threads
  THREADS="$2"
  shift
  ;;
  *)
  
  ;;
esac
shift
done

PATH="/usr/local/bin/anaconda3/bin:$PATH"

echo RUNNING FASTQ CONVERSION ON ALL .SRA FILES IN "${INDIR}" WITH "${THREADS}" THREADS
echo OUTPUT DIRECTORY = "${OUTDIR}"
INFILES="$INDIR/SRR*.sra"
for f in $INFILES
do
  echo
  echo CONVERTING $f
  echo
  SRA_ID=${f%%.sra}
  parallel-fastq-dump --sra-id $f -t $THREADS -O $OUTDIR --tmpdir $INDIR --split-files
done
