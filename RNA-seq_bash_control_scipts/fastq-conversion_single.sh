#!/bin/bash

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
  -i|--indir) # directory for sra files
  INDIR="$2"
  shift
  ;;
  -O|--outdir) # directory to write fastq files
  OUTDIR="$2"
  shift
  ;;
  --suffix) # suffix of files to convert
  SUFFIX="$2"
  shift
  ;;
  *)

  ;;
esac
shift
done

# get a list of sra files from INDIR with SUFFIX
cd $INDIR
INFILES=*$SUFFIX

# loop over files to trim
for f in $INFILES
do
  echo
  echo Converting $f
  echo
  OUTPATH="$OUTDIR/"
  /mnt/disks/disk-2/Tools/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump.2.8.2 $f -O $OUTPATH
done
