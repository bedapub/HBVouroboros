#/bin/bash

# Extract sequences from BAM files of paired reads into two FASTQ files, and run Trinity
# zhangj83

function help {
    echo "Extract sequence from one or more BAM files of paired reads, write them into two FASTQ files, and run Trinity"
    echo "Usage: $0 [-o|--outdir OUTDIR] [-c|--cpu CPU] [-h|--help] bamfiles..."
    echo
    echo "Mandatory parameters:"
    echo "  bamfiles: BAM files, wild cards are accepted"
    echo
    echo "Optional parameters:"
    echo "  -h|--help: print this message and quit"
    echo "  -c|--cpu: number of threads (CPUs) used"
    echo "  -o|--outdir: Set output directory. Default output directory: TrinityFromBams_outdir"
    echo
    echo "Contact: Jitao David Zhang <jitao_david.zhang@roche.com>, Tel 86251"
}

## CONSTANTS
## parse parameters
OUTDIR="TrinityFromBams_outdir"
cpu=1
declare -a bamfiles
while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in 
      -o|--outdir)
	  OUTDIR="$2"
	  shift ## past argument
	  shift ## past value
	  ;;
      -c|--cpu)
	  cpu="$2"
	  shift
	  shift
	  ;;
      -h|--help)
	  help
	  exit 0
	  ;;
      *)
	  bamfiles+=("$1")
	  shift
	  ;;
  esac
done

if [ ${#bamfiles[@]} -eq 0 ]; then
    help
    exit 0
fi

## modules
ml load SAMtools
ml load Bowtie

## Implement logic
mkdir -p "$OUTDIR"

fq1="$OUTDIR"/mergedReads-1.fq.gz
fq2="$OUTDIR"/mergedReads-2.fq.gz

rm -rf ${fq1} ${fq2}

mergedBam="$OUTDIR"/merged.bam

samtools merge -@ "$cpu" "$mergedBam" ${bamfiles[@]}
samtools fastq -@ "$cpu" -c 9 "$mergedBam" -1 "$fq1" -2 "$fq2"

## run trinity
/pstore/home/zhangj83/apps/trinity/trinityrnaseq-2.2.0/Trinity --seqType fq --max_memory 16G --left "$fq1" --right "$fq2" --output "$OUTDIR"/trinity_output_dir

## program exist
