#!/bin/bash

## This shows how we estimate the local mutation rate in the human genome in
## 100-kb windows

## We do this by downloading multiple alignments from the UCSC Genome Browser
## and running the phyloP software in every window.

## It requires phyloP and tree_doctor utilities which are part of the PHAST software
## (http://compgen.cshl.edu/phast)

## It also requires bedops (https://bedops.readthedocs.io/en/latest)
## and the tool "mafRanges" which is available as part of the UCSC genome browser
## software:
## (http://hgdownload.soe.ucsc.edu/admin/exe/<OS_directory>)


## First download global substitution rates and a phylogenetic tree fit to
## the 100-way human alignments
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.mod
seqs="hg19,panTro4,gorGor3,ponAbe2,nomLeu3"
tree_doctor --prune-all-but $seqs hg19.100way.phastCons.mod > primates.mod

## we are going to mask out regions within $flank of conserved elements
## since we are interested in neutral rates
flank=100

## These files contain predicted conserved elements various subsets of the alignments:
for subset in 46way 46wayPlacental 46wayPrimates 100way; do
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/phastConsElements$subset.txt.gz
done


for chrom in `seq 1 22`; do

    ## download and unzip the MAF file
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/multiz100way/maf/chr$chrom.maf.gz
    gunzip chr$chrom.maf.gz
    mafFile=chr$chrom.maf

    ## get coordinate range of alignment
    mafRanges $mafFile hg19 /dev/stdout |
       awk -v OFS="\t" 'NR==1 {firstCoord=$2};
                        END{print $1,firstCoord,$3}' > mafRange$chrom.bed

    ## Make windows of size 100000; will estimate a substitution rate in every window
    bedops --chop 100000 --stagger 10000 mafRange$chrom.bed > overlapping_windows$chrom.bed
    bedops --chop 10000 mafRange$chrom.bed > windows${chrom}.bed

    lastCoord=`tail -n 1 mafRange$chrom.bed | cut -f 3`
    ## get phastCons elements on this chromosome with flanking regions on either end

    ## Get the regions of the MAF to keep. We want to mask all the phastConsElements
    ## with flanking regions on either side.
    ## This loop combines all those regions from the phastConsElements tables,
    ## then takes the difference between the whole mafRegion and these.
    ( for f in phastConsElements*.txt.gz; do
        zcat $f |
          awk -v OFS="\t" -v chrom=chr$chrom -v lastCoord=$lastCoord -v flank=$flank '
             $2 == chrom {
                   chromStart=$3-flank;
                   chromEnd=$4+flank;
                   if (chromStart < 0) {chromStart=0};
                   if (chromEnd > lastCoord) {chromEnd=lastCoord};
                   print chrom,chromStart,chromEnd}'
      done ) | sort-bed - |
        bedops -d mafRange$chrom.bed -  > keepRegions$chrom.bed

    ## Now filter the maf file to keep only those regions, only keeping
    ## some closely related primates (human, chimp, gorilla, orangutan, gibbon)
    maf_parse --features keepRegions$chrom.bed \
	      --seqs "hg19,panTro4,gorGor3,ponAbe2,nomLeu3" $mafFile > chr$chrom.filtered.maf


    ## run phyloP (this takes hours for some chromosomes)
    phyloP -i MAF --method LRT --features overlapping_windows$chrom.bed --mode CONACC \
	   hg19.100way.phastCons.mod chr$chrom.filtered.maf > phyloP_${chrom}.txt

    ## convert phyloP scores to bed file with columns chrom, chromStart, chromEnd, name (not used), score
    ## Then use bedmap to take average "score" in each window
    ## The score is the "scale" parameter estimated by phyloP which gives the substitution rate relative
    ## to neutrality
    ## In some cases, there is no alignment data in the window that has not been masked.
    ## Usually these are due to big gaps in the middle of an alignment (often near centromere)
    ## In these cases set relative rate to 1.0. This occurs in 1.6% of windows.
    cut -f 1-5 phyloP_${chrom}.txt | grep -v "^#" |
	bedmap --echo --mean --delim "\t" windows${chrom}.bed - |
	sed 's/NAN/1.0/g' > relSubstRates_$chrom.bed

done

## Once all chroms are done, combine into one file and scale so that average value
## Is mean mutation rate

## First get current average value; can just average 5th column since all
## windows are equal in size
meanVal=`cat relSubstRates_*.bed |
   awk '$0 !~ /^#/ {s+=$4; total++}; END{print s/total}'`

## Now create a bed file whose fourth column is scaled rate
cat relSubstRates_*.bed | sed 's/chr//g' | sort-bed - |
  awk -v newMean=1.45e-8 -v oldMean=$meanVal -v OFS="\t" '
     $0 !~ /^#/ {print $1,$2,$3,$4*newMean/oldMean}' |
      bgzip > subst_rate_autosome.bed.gz
