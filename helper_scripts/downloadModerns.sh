#!/bin/bash

## Download local VCFs for 14 modern humans


if [[ $# -ne 2 ]]; then
    echo "Usage: bash downloadModerns.sh <outdir> <ind> <chr:start-end>"
    exit 1
fi

outdir=$1
region=$2

baseURL=http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf


if [[ ! -e ModernHumans/README_Bteam ]]; then
    mkdir -p ModernHumans
    cd ModernHumans
    wget $baseURL/README_Bteam
    cd ..
fi


inds=`grep -o "SS[0-9]*" ModernHumans/README_Bteam`
chrom=`echo $region | cut -f 1 -d ':'`
mkdir -p $outdir

for ind in $inds; do
    ## Rename each ind by population
    indName=`grep $ind ModernHumans/README_Bteam |
       cut -f 2`
    ## There are two Australians so give them a number
    if [[ $indName == "Australian" ]]; then
       if [[ $ind == "SS6004477" ]]; then
          indName="Australian.1"
       else
           indName="Australian.2"
       fi
    fi
    tabix -h $baseURL/$ind.hg19_1000g.$chrom.mod.vcf.gz $region |
       awk -v ind=$indName -v OFS="\t" '
          $0 ~ /^#CHROM/ && NF==10 {$10=ind};
          {print $0}' | bgzip > $outdir/$indName.vcf.gz
    tabix -p vcf $outdir/$indName.vcf.gz
done

