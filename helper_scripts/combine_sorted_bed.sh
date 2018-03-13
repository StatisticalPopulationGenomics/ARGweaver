#!/bin/bash

# take two sorted bedGraph-type files and return a merged file with both sets of
# scores. Each bed file should have only one row for any genomic region, but can
# have an arbitrary number of columns (though the same throughout file).
# So if bedfile1 is chrom,chromStart,chromEnd,score1,score2,score3
# and bedfile2 is   chrom,chromStart,chromEnd,score4,score5
# Then resulting bedfile will have all 5 scores on each row
# if bedfiles are gzipped then will unzip them to a temporary directory.
# Bedfiles cannot be streamed into this script, as they are read multiple times

bedfile1=$1
bedfile2=$2

tmpdir=`mktemp -d`
isZipped=`file $bedfile1 | grep gzip`
if [[ -n $isZipped ]]; then
  gunzip -c $bedfile1 > $tmpdir/file1.bed
  bedfile1=$tmpdir/file1.bed
fi
isZipped=`file $bedfile2 | grep gzip`
if [[ -n $isZipped ]]; then
  gunzip -c $bedfile2 > $tmpdir/file2.bed
  bedfile2=$tmpdir/file2.bed
fi
numcol1=`awk -F "\t" '$0 !~ /^#/ {print NF; exit}' $bedfile1`
numcol2=`awk -F "\t" '$0 !~ /^#/ {print NF; exit}' $bedfile2`


bedops --partition --header $bedfile1 $bedfile2 |
  bedmap --header --delim "\t" --echo --echo-map - $bedfile1 |
  bedmap --header --delim "\t" --echo --echo-map - $bedfile2 |
  cut -f 1-3,7-$((3+$numcol1)),$((7+$numcol1))-$((3+$numcol1+$numcol2))

