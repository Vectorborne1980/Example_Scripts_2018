#!/bin/bash

function usage {
  echo "Usage: $0 [--indir] [--outdir] [--suffix] [--threads] [--ilmnclip] [--headcrop] [--slidingwindow]"
  echo
  echo "Runs Trimmomatic on all .fastq files in a specified input directory and writes trimmed .fastq files to an output directory."
  echo
  echo "  -i/--indir        input directory of fastq files"
  echo "  -o/--outdir       output directory for trimmed fastq files"
  echo "  -t/--threads      number of threads"
  echo "  --ilmnclip        argument for ILLUMINACLIP. This includes a path to the desired adapter file, followed immediately by :seedMismatches:palindromeClipThreshold:simpleClipThreshold (:2:30:10)"
  echo "  --headcrop        argument for HEADCROP, usually about 14 biased bases at the start of illumina reads"
  echo "  --slidingwindow   argument for SLIDINGWINDOW. windowSize:requiredQuality e.g. 4:30"
  echo "  See Trimmomatic documentation for more help: www.usadellab.org/cms/?page=trimmomatic"
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
  -t|--threads) # number of threads
  THREADS="$2"
  shift
  ;;
  -i|--indir) # directory for fastq files
  INDIR="$2"
  shift
  ;;
  -o|--outdir) # directory to write trimmed files
  OUTDIR="$2"
  shift
  ;;
  --ilmnclip) # ILLUMINACLIP argument.  /path/adapter.fa:seedmismatches:palindromeclipthreshold:simpleclipthreshold (e.g. /path:2:30:10).  Still must specify all parameters even in single end mode. 
  ILMNCLIP="$2"
  shift
  ;;
  --headcrop) # HEADCROP parameter.  usually 14 for illumina
  HEADCROP="$2"
  shift
  ;;
  --slidingwindow) # SLIDINGWINDOW argument.  windowsize:requiredquality e.g. 4:30
  SLIDINGWINDOW="$2"
  shift
  ;;
  *)
  
  ;;
esac
shift
done

# get a list of -basein files from INDIR with SUFFIX
cd $INDIR
INFILES=SRR*.fastq
SUFFIX=.fastq

# loop over files to trim
for f in $INFILES
do
  echo
  echo TRIMMING $f
  echo
  OUTFILE=${f%%$SUFFIX}"_trimmed.fastq"
  OUTPATH="$OUTDIR/$OUTFILE"
  java -jar '/mnt/disks/disk-2/Tools/Trimmomatic-0.36/trimmomatic-0.36.jar' SE -threads $THREADS $f $OUTPATH ILLUMINACLIP:$ILMNCLIP HEADCROP:$HEADCROP SLIDINGWINDOW:$SLIDINGWINDOW
done

