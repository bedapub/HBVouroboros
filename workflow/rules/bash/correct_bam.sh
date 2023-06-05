#!/bin/bash

while getopts f:s:o: arg; do

        case "${arg}" in
                f) FastaInput=${OPTARG}
                        ;;
                s) sortBamInput=${OPTARG}
                         ;;
                o) outputFile=${OPTARG} 
			;;
        esac
done


line=$(awk '{if(NR==1) print $0}' ${FastaInput} | cut -d"=" -f2 | cut -d" " -f1)
echo $line

samtools view -h ${sortBamInput} | awk -v Line="$line" -v FS="\t" -v OFS="\t" '{{if ($0 ~ "^[^@]" && $4>Line) {$4 =$4-Line;}} print $0 }' | samtools sort | samtools view -b > ${outputFile};
