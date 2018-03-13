#!/bin/bash

# This script documents how the recombination map file was created.
# It downloads recombination maps from Hinch et al 2011 based on
# African American genomes. It then smooths it over 5kb windows and
# converts it to a 4-column sorted bed file.

wget http://www.well.ox.ac.uk/~anjali/AAmap/maps_b37.tar.gz
tar -zxf maps_b37.tar.gz

## This file contains several recombination maps in hg19;
## I will use the African American (AA) one
## as it appears to fall in the middle of the range and as an African-based
## map should be more appropriate for an analysis of deep human ancestry than
## a European-based map.
## The physical position is in the 1st column
## The AA genetic map position (in cM) is in the 6th column

## get average rate per generation per base-pair by computing
## (change in genetic position)/(change in physical position)*0.01
## (since centimorgans represent a 1% probability of recombination)
##
## This command creates a bed file with average rate between every two points
## in the map
rm -f recomb_map.AArates.bed
for chrom in `seq 1 22 | sort`; do
   awk -v OFS="\t" -v chrom=$chrom '
      NR > 2 {print chrom,prevPhyPos,$1,"AArate",($6-prevGenPos)/($1-prevPhyPos)*0.01};
     {prevGenPos=$6;
      prevPhyPos=$1}' maps_b37/maps_chr.$chrom >> recomb_map.AArates.bed
done

## Now smooth map in desired region.
# First make windows of size 5kb
bedops -m recomb_map.AArates.bed |
    bedops -w 5000 - > windows.bed

## Now get mean within each window. Impose a minimum rate of 2e-10
bedmap --echo --wmean --sci --delim "\t" windows.bed recomb_map.AArates.bed |
    awk -v OFS="\t" '$4 < 2e-10 {$4=2e-10}; {print $0}' |
    bash merge_bed.sh |
    bgzip > recomb_rate_autosome.bed.gz

## clean up
rm -rf windows.bed maps_b37 maps_b37.tar.gz recomb_map.AArates.bed
