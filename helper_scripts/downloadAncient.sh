## This script shows how the ancient sample VCF files are obtained
## The ancient samples are described in these papers:
## Altai Neandertal (from Prufer et al 2014, "The complete genome sequence
##     of a Neanderthal from the Altai Mountains, Nature)
## Denisovan (from Meyer et al 2012, "A high-coverage genome sequence
##     from an Archaic Denisovan Individual, Science)
## Vindija Neandertal (from Prufer et al 2017 "A high-coverage Neandertal genome
##     from Vindija Cave in Croatia", Science)

## The Altai and Denisovan Neandertals SNP calls were recomputed for the
## most recent Vindija paper so that all three sets were done consistently,
## so we will use that set of calls.

if [[ $# -ne 3 ]]; then
    echo "Usage: bash downloadAncient.sh <CHR:START-END> <ind> <OUT.vcf.gz>"
    exit 1
fi

region=$1
ind=$2
outFile=$3

# extract chromosome name from region
chrom=`echo $region | cut -f 1 -d ':'`

ancientURL=http://cdna.eva.mpg.de/neandertal

if [[ $ind == "Altai" ]]; then
    file=Vindija/VCF/Altai/chr${chrom}_mq25_mapab100.vcf.gz
elif [[ $ind == "Vindija" ]]; then
    file=Vindija/VCF/Vindija33.19/chr${chrom}_mq25_mapab100.vcf.gz
elif [[ $ind == "Denisova" ]]; then
    file=Vindija/VCF/Denisova/chr${chrom}_mq25_mapab100.vcf.gz
fi

## Do the download in a different directory for each individual
## because the file names are the same (but in different directories)
## for each ind on the server.  tabix downloads index files and
## this prevents it using the same index file for different individuals
downloadDir=${ind}_index_files
mkdir -p $downloadDir
cd $downloadDir
tmpout=`mktemp --tmpdir=.`

## Download the file and rename the individual to $ind
tabix -f -h $ancientURL/$file $region |
 awk -v ind=$ind -v OFS="\t" '
    $0 ~ /^#CHROM/ && NF==10 {$10=ind};
           {print $0}' | bgzip > $tmpout
cd ..
mv $downloadDir/$tmpout $outFile

## create index file for small VCF that was dowloaded
tabix -p vcf $outFile
