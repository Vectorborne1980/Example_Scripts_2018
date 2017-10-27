#!/bin/bash

function usage {
  echo "Usage: $0 [--indir] [--outdir] [--stranded] [-M] [--threads]"
  echo
  echo "Runs featureCounts on sorted bam files and writes countsall.txt and countsall.txt.summary reports to an output directory. Requires both ends of paired end reads to map (-B option in featureCounts)"
  echo
  echo "  -i/--indir    input directory of sorted bam files"
  echo "  -o/--outdir   output directory for counts files"
  echo "  -s/--stranded 0 for unstranded, 1 for stranded, 2 for reverse stranded"
  echo "  -M            1 to count multimapping reads, 0 to discard (default)"
  echo "  -t/--threads  number of threads"
  echo "  Version control: Nick Geraci July 27, 2017"
  echo
  exit 1
}

if [ "$1" == "-h" ]
then
  usage
fi

# setting up mutlimapping reads option.  False by default
M=0

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
  -i|--indir) # input directory of sorted bam files
  INDIR="$2"
  shift
  ;;
  -o|--outdir) # output directory
  OUTDIR="$2"
  shift
  ;;
  -s|--stranded) # 0 for unstranded, 1 for stranded, 2 for reverse stranded
  STRANDED="$2"
  shift
  ;;
  -t|--threads)
  THREADS="$2"
  shift
  ;;
  -M) # count multimapping reads
  M="$2"
  shift
  ;;
  *)
  
  ;;
esac
shift
done

INFILES="$INDIR/*_sorted"
OUTPATH="$OUTDIR/countsall.txt"

echo
echo COUNTING FEATURES IN SORTED BAM FILES IN $INDIR
echo
if [[ $M = 1 ]]
then
/mnt/disks/disk-2/Tools/subread-1.5.3-Linux-x86_64/bin/featureCounts -p -s $STRANDED -T $THREADS -M -B -a /mnt/disks/disk-2/Tools/Reference_Genomes/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf -F GTF -t exon -g gene_id -o $OUTPATH $INFILES
else
/mnt/disks/disk-2/Tools/subread-1.5.3-Linux-x86_64/bin/featureCounts -p -s $STRANDED -T $THREADS -B -a /mnt/disks/disk-2/Tools/Reference_Genomes/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf -F GTF -t exon -g gene_id -o $OUTPATH $INFILES
fi
