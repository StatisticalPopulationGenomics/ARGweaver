#!/bin/bash

## This script walks through an ARGweaver analysis for a genomic interval,
## including downloading data, creating masks and parameter files,
## running ARGweaver, and doing some post-analysis.

## It is meant to be run in a bash shell on linux, though should also work
## on a Mac.
## For Windows users, it is suggested to run through a linux Virtual Machine
## (for example using Oracle VirtualBox)

## Besides bash, the following software will need to be installed and in your PATH:
## ARGweaver (https://github.com/CshlSiepelLab/argweaver.git)
## samtools (www.htslib.org)
## bedops (https://bedops.readthedocs.io/en/latest)
## phast (http://compgen.cshl.edu/phast/): Note this is only required for computing
##    regional substitution rates- for the human genome (hg19 assembly) these
##    are provided for autosomal chromosomes
## bigWigToBedGraph (http://hgdownload.soe.ucsc.edu/admin/exe/<OS_directory>)

##############################################################
########### PART 0: set variables/directory paths ############
##############################################################

## Set this to argweaver installation directory
ARGWEAVER_HOME=~/argweaver

export PATH=$ARGWEAVER_HOME/bin:$PATH
export PYTHONPATH=$ARGWEAVER_HOME:$PYTHONPATH

region=1:159113000-159238000
regionName=duffy

## create a directory for the files
mkdir -p $regionName

# maxtime is the maximum time (in generations) in ARGweaver's set of
# discrete time points
maxtime=200e3

## delta determines the spacing between discrete time points- 0.01 is the
## default. Smaller numbers produce closer spacing at recent time intervals
delta=0.01


## parse region to get chromosome and region bounds
chrom=`echo $region | cut -f 1 -d ':'`
chromStart=`echo $region | awk -F "[:-]" '{print $2}'`
chromEnd=`echo $region | awk -F "[:-]" '{print $3}'`


#######################################################################
################ PART 1: Download VCF Files ###########################
#######################################################################

for ind in Altai Denisova Vindija; do
    if [[ ! -e $regionName/$ind.vcf.gz ]]; then
        bash helper_scripts/downloadAncient.sh $region $ind $regionName/$ind.vcf.gz
    fi
done

bash helper_scripts/downloadModerns.sh $regionName $region


## Assign inds to population groups and create files for each group
afrInds="Mandenka Mbuti San Yoruba Dinka"
eurAsiInds="French Sardinian Karitiana Han Papuan Australian.1 Australian.2 Mixe Dai"
ancientInds="Altai Vindija Denisova"
inds="$afrInds $eurAsiInds $ancientInds"

( for ind in $afrInds; do echo $ind; done ) > $regionName/african_inds.txt
( for ind in $eurAsiInds; do echo $ind; done ) > $regionName/eurasian_inds.txt
grep -v San $regionName/african_inds.txt > $regionName/african_no_san_inds.txt
(for ind in $afrInds $eurAsiInds; do echo $ind; done ) > $regionName/human_inds.txt
( for ind in $ancientInds; do echo $ind; done ) > $regionName/ancient_inds.txt
( for ind in $inds; do echo $ind; done ) > $regionName/inds.txt

## make list of VCF files
ls $regionName/*.vcf.gz > $regionName/vcf_files.txt


#######################################################################
################ PART 2: Create mask files ############################
#######################################################################

## mask out regions that have been flagged for poor mapability/alignment
if [[ ! -e mask_files/filter.bed.gz ]]; then
    mkdir -p mask_files
    cd mask_files

    ## download blacklist tables from UCSC
    blackListTables="wgEncodeDacMapabilityConsensusExcludable wgEncodeDukeMapabilityRegionsExcludable"
    ucscDir=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database
    mkdir -p blackList_tmp
    for t in $blackListTables; do
        wget $ucscDir/$t.txt.gz
        gunzip -c $t.txt.gz | awk -v OFS="\t" '{print $2,$3,$4}' |
            sed 's/^chr//g' | sort-bed - > blackList_tmp/$t.bed
    done
    bedops -m blackList_tmp/*.bed > blackList.bed
    rm -rf blackList_tmp

    # This requires bigWigToBedGraph utility, available on UCSC's webpage:
    # http://hgdownload.soe.ucsc.edu/admin/exe/<OS_directory>
    if [[ ! -e mapability75_1.bed.gz ]]; then
      ucscDir=http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi
      bigWigToBedGraph $ucscDir/wgEncodeCrgMapabilityAlign75mer.bw /dev/stdout |
         awk '$4 != 1' | sed 's/chr//g' | sort-bed - | bedops -m - \
          | bgzip > mapability75_1.bed.gz
    fi
    gunzip mapability75_1.bed.gz
    bedops -m blackList.bed mapability75_1.bed | bgzip > filter.bed.gz
    rm -f blackList.bed
    bgzip mapability75_1.bed
    cd ..
fi


## The VCF files we have downloaded contains a row for every site with data,
## (this is different from most VCF files which contain only confident
## variants).
## Therefore any missing rows should be treated as missing data;
## create a filter for each individual containing sites absent from VCF
echo -e "$chrom\t$(($chromStart-1))\t$chromEnd" > $regionName/region.bed
rm -f $regionName/ind_mask_files.txt
for ind in $inds; do
     gunzip -c $regionName/$ind.vcf.gz |
        awk '$0 !~ /^#/ {print $1,$2-1,$2}' |
        bedops -d $regionName/region.bed - |
        bgzip > $regionName/$ind.mask.bed.gz
     echo -e "$ind\t$regionName/$ind.mask.bed.gz" >> $regionName/ind_mask_files.txt
done


#######################################################################
################ PART 3: Create mutation rate file ####################
#######################################################################

## Make substitution rate map. This is slow, so only do once per
## chromosome. The maps have been pre-computed for all chromosomes of hg19
## and are in the directory subst_rates

## This requires the software PHAST, available at:
## http://compgen.cshl.edu/phast
## which provides the utilities phyloP and tree_doctor

## subst_rate file is created by script helper_scripts/subst_rate.sh
## It takes awhile so just get pre-computed rates

if [[ ! -e rate_files/subst_rate_autosome.bed.gz ]]; then
    echo "Subst rate file not found" >> /dev/stderr
    exit 1
fi



#######################################################################
################ PART 4: Create recombination rate file ###############
#######################################################################

## Get African American map from Hinch et al:
# "The landscape of recombination in African Americans" Nature, 2011

## This is done in separate script helper_scripts/recomb_rate.sh
## The precomputed rates are in github
if [[ ! -e rate_files/recomb_rate_autosome.bed.gz ]]; then
    echo "Recomb rate file not found" >> /dev/stderr
    exit 1
fi


#######################################################################
################ PART 5: Choose population size #######################
#######################################################################

# In this example, samples are coming from Africa, Europe, Asia, as well
# as ancient populations. It is impossible to choose a demographic history
# that will fit this range of individuals. So we will use the default of
# 10000, which is a very rough estimate of human effective population size.

# When individuals come from a single population, a more realistic
# demographic model may produce better ARGs. A 2 column file with rows
# <time1> <popsize1>
# <time2> <popsize2>
# ...
# indicates that the population has diploid popsize1 until time1
# (where time1 is the number of generations before present)
# then the population size changes to popsize2 until time2, etc
# Then this file can be passed to arg-sample with the option --popsize-file


#######################################################################
################ PART 6: Set ancient sample ages ######################
#######################################################################

## Also create sample_ages.txt, which gives the age (in generations)
## of the ancient genomes.
## Use the dates computed in Prufer et al 2017:
## Altai: 122000 years
## Denisvova: 72000 years
## Vindija: 52000 years
## assume 29 year generations
if [[ ! -e sample_ages.txt ]]; then
( echo -e "Altai\t$((122000/29))"
  echo -e "Denisova\t$((72000/29))"
  echo -e "Vindija\t$((52000/29))"
) > sample_ages.txt
fi



#######################################################################
################ PART 7: Run ARGweaver ################################
#######################################################################

## This is another potential filter but it's pretty heavy handed
#bigWigToBedGraph -chrom=chr$chrom -start=$startCoord -end=$endCoord \
#   http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign75mer.bw \
#   /dev/stdout  | awk '$4 != 1' | sed 's/chr//g' | bedops -m | bgzip \
# > uniqueness.bed.gz

outdir=$regionName/args
mkdir -p $outdir

arg-sample --vcf-files $regionName/vcf_files.txt \
   --region $region \
   --vcf-min-qual 30 \
   --subsites $regionName/inds.txt \
   --maskmap mask_files/filter.bed.gz \
   --mask-cluster 2,5 \
   --ind-maskmap $regionName/ind_mask_files.txt \
   --age-file sample_ages.txt \
   --mutmap rate_files/subst_rate_autosome.bed.gz \
   --recombmap rate_files/recomb_rate_autosome.bed.gz \
   -c 10 \
   --maxtime $maxtime \
   --delta $delta \
   -n 2000 \
   -o $outdir/out $overwrite


#######################################################################
################ PART 8: Post analysis ################################
#######################################################################

## This section creates the "layout" file for creating
## leaf trace plots. The first two commands just determine input and
## output file names, based on the last SMC file written by arg-sample.
## The "arg-layout" command creates the layout file.
lastSmc=`ls -t $outdir/*.smc.gz | head -n1`
layout=`echo $lastSmc | sed 's/smc.gz$/layout.gz/'`
echo $layout $lastSmc
arg-layout $layout $lastSmc

## Set burnin
burnin=600

## This converts all the ARGss into a single indexed bed file
## (which is the format required by arg-summarize program)
smc2bed-all $outdir/out

## Compute median stats and confidence intervals for each population
## for TMRCA, branchlen, RTH, popsize, pi, tmrca-half
for pop in all human eurasian african; do
    if [[ $pop == "all" ]]; then
        subsetArg="--subset-inds $regionName/inds.txt"
    else
        subsetArg="--subset-inds $regionName/${pop}_inds.txt"
    fi
    for stat in tmrca branchlen rth popsize pi tmrca-half; do
        arg-summarize -a $outdir/out.bed.gz --log-file $outdir/out.log \
           --$stat --quantile 0.05,0.5,0.95 --burnin $burnin $subsetArg |
           bgzip > $outdir/out.$stat.$pop.bed.gz
    done
    arg-summarize -a $outdir/out.bed.gz --log-file $outdir/out.log $subsetArg \
       --branchlen --node-dist-all --burnin $burnin --quantile 0.5 |
        bgzip > $outdir/out.node_dists.$pop.bed.gz
done

## This script simply combines results across different sub-populations into
## a single file
for stat in tmrca branchlen rth popsize pi tmrca-half max-coal-rate; do
    echo $stat
    for pop in eurasian african african_no_san; do
	bash helper_scripts/combine_sorted_bed.sh $outdir/out.$stat.$pop.bed.gz \
	       $outdir/out.$stat.human.bed.gz |
           bgzip > $outdir/out.$stat.$pop.human.bed.gz
done
done



## This command computes median distances between every haplotype pair at every
## local tree.
## (It is commented out because it creates very large output)
#arg-summarize -a $outdir/out.bed.gz --log-file $outdir/out.log --subset-inds $regionName/inds.txt \
#   --branchlen --node-dist-all --burnin $burnin --quantile 0.5 |
#  gzip > $outdir/out.node_dists.bed.gz


## This command gets minimum coalescence time between each human individual
## and each ancient individual
humanInds=`cat $regionName/human_inds.txt`
for ind in $humanInds; do
    arg-summarize -a $outdir/out.bed.gz --log-file $outdir/out.log \
       --min-coal-time "$ind,Altai;$ind,Vindija;$ind,Denisova" \
       --quantile 0.05,0.5,0.95 --burnin $burnin |
      bgzip > $outdir/out.$ind.ancient.bed.gz
done


## This computes the average allele age for all SNPs at each rep
## The -B argument tells it to first prune ancient individuals
allele_age -s $burnin -B $regionName/human_inds.txt \
	   -t $regionName/tmp $outdir/out


