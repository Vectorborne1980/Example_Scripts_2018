#!/bin/bash

function usage {
  echo "Usage: $0 [--indir] [--outdir] [--threads]"
  echo
  echo "Runs FastQC on all .fastq files in a specified input directory and writes .html reports to an output directory. Also compiles a MultiQC report in the output directory. Can run as many files in parallel as there are threads."
  echo
  echo "  -i/--indir    input directory of fastq files"
  echo "  -o/--outdir   output directory for html reports"
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
  -i|--indir) # directory to read fastq files from
  INDIR="$2"
  shift
  ;;
  -o|--outdir) # directory for output html reports
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

echo RUNNING FASTQC ON ALL .FASTQ FILES IN "${INDIR}" WITH "${THREADS}" THREADS
echo OUTPUT DIRECTORY = "${OUTDIR}"
INFILES="$INDIR/*.fastq"
'/mnt/disks/disk-2/Tools/FastQC/fastqc' $INFILES --threads $THREADS --outdir $OUTDIR
cd $OUTDIR
multiqc .
