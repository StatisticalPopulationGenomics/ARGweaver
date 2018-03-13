require("argweaver")


## Read in an output file from arg-summarize
## assume it has columns named <stat>_quantile_0.500,
## <stat>_quantile_0.050, and <stat>_quantile_0.950
## plots bigWig style plot
makeBwPlot <- function(file, stat, scale=1, ylim=NULL, doBounds=TRUE, ...) {
    if (stat=="tmrca-half") stat <- "tmrca_half"
    x <- readArgSummary(file)
    medianName <- paste(stat, "_quantile_0.500", sep="")
    if (!doBounds) {
            if (is.null(ylim)) ylim <- c(0, max(x[,medianName]*scale))
            plotBw(x, x[,medianName]*scale, ylim=ylim, ...)
    } else {
        lowName <- paste(stat, "_quantile_0.050", sep="")
        highName <- paste(stat, "_quantile_0.950", sep="")
        if (is.null(ylim)) ylim <- c(0, max(x[,highName]*scale))
        plotBw(x, x[,medianName]*scale, y0=x[,lowName]*scale, y1=x[,highName]*scale, ylim=ylim, ...)
    }
}                



region <- commandArgs(trailingOnly=TRUE)[1]
if (is.na(region))
    region <- "duffy"

## change to region's directory; there should be a sub-directory named "args"
## with argWeaver results
setwd(region)

## Set up list assigning color to each individual based on population
popColors <- list()
popColors[["all"]] <- "black"
popColors[["eurasian"]] <- "coral3"
popColors[["african"]] <- "cyan4"
popColors[["african_no_san"]] <- "blue"
popColors[["ancient"]] <- "black"
popColors[["human"]] <- "orange"
pops <- names(popColors)
indColors <- list()
for (pop in pops) {
    file <- paste0(pop, "_inds.txt")
    if (file.exists(file)) {
        inds <- scan(file, what=character())
        if (pop != "human") indColors[inds] <- popColors[[pop]]
    }
}


## Get region boundaries and genes within region
boundaries <- scan("region.bed")
chrom <- boundaries[1]
regionStart <- boundaries[2]+1 ## make 1-based
regionEnd <- boundaries[3]

## Change to output arg directory
setwd("args")

## Plot basic statistics
par(cex.lab=1.4, cex.axis=1.2)
plotTraces("out.stats", stats=c("prior", "likelihood", "joint", "recombs", "arglen", "noncompats"))


## Or, look at individual statistics, such as "joint"
par(mfrow=c(1,1))
x <- readStatsFile("out.stats")
plotStats(x, "joint")


## Plot autocorrelation after burnin
acf(x[x$iter > 1000,"joint"])

## Look at leafTrace
file <- list.files(".", pattern=".*.layout.gz$")
plotLeafTrace(file, col=indColors)



## make plots of each statistic for each population, with confidence intervals
par(mgp=c(2.5,1,0), mfrow=c(1,1), mar=c(4,4,2,1), mfrow=c(4,1))
xlim <- c(159.13e6, 159.22e6)
for (stat in c("pi", "tmrca", "rth", "popsize")) {
    add <- FALSE
    if (stat=="rth") {
        scale <- 1
        ylab <- "RTH"
        ylim <- c(0,1)
    } else if (stat == "popsize") {
        scale <- 1
        ylab <- "Popsize"
    } else {
        scale <- 29/1000
        ylab <- paste(stat, "(kya)")
    }
    ylim <- NULL
    for (pop in c("eurasian", "african", "african_no_san")) {
        file <- paste0("out.", stat, ".", pop, ".bed.gz")
        if (file.exists(file)) {
            makeBwPlot(file, stat, add=add, col=popColors[[pop]], ylab=ylab, scale=scale, ylim=ylim, doBounds=TRUE,
                       xlim=xlim)
            add <- TRUE
        }
    }
}



## Plot a single tree
plotTreesFromBed("out.bed.gz", chrom=chrom, start=159173860, end=159173860, leafCol=indColors,
                 regionLine=5, timeScale=29/1000, ylab="Kya")

## Here is a command to browse trees every 5kb
par(ask=TRUE)  ## This is good for interactive mode; don't show next plot until ENTER is pressed
plotTreesFromBed("out.bed.gz", chrom=chrom, start=xlim[1], end=xlim[2], leafCol=indColors, regionLine=5,
                 interval=5000, timeScale=29/1000, ylab="Kya", ylim=c(0, 2000))


## This creates a plot of sites in introgressed region
x<- readSites("out.2000.sites.gz")
indOrder=c("Denisova", "Altai", "Vindija", "Han", "Dai", "Yoruba", "San")
indOrder <- sprintf("%s_%i", sapply(indOrder, rep, 2), 2:1)
plotSites(x=x, start=159196000, end=159209000, stretch=20,
          indOrder=indOrder, textColor=indColors, useCoords=FALSE,
          xlab="Variant Site Index")



### Allele age analysis
file <- "out.allele_age.bed.gz"
x <- readArgSummary(file)
sum(x$inf_sites == 0) / nrow(x)   ### This is fraction of sites requiring multiple mutations

## Remove sites that require multiple mutations
x <- x[x$inf_sites==1,]

## Use the coordinate as a SNP identifier
snpID <- as.character(x$chromEnd)

## For each unique SNP, get the most frequent derived allele
derAllele <- tapply(as.character(x$derAllele), snpID, function(x) {names(which.max(table(x)))})

## This produces a logical vector for each snp/sample saying whether the
## inferred derived allele is different from the consensus derived allele
flipped <- derAllele[snpID] != x$derAllele

## Look at how many are "flipped" as a function of allele frequency
freqs <- sort(unique(x$derFreq))
unFlippedCounts <- table(x[!flipped,"derFreq"])[freqs]
totalCounts <- table(x[,"derFreq"])[freqs]

## Plot the fraction of
plot(freqs, as.numeric(unFlippedCounts/totalCounts), ylim=c(0,1),
     ylab = "Fraction allele flipped", xlab="Derived frequency")

## Now look only at allele age estimates across reps with the same derived allele
x <- x[!flipped,]
snpID <- snpID[!flipped]


x$freq <- x$derFreq/(x$derFreq + x$ancFreq)
medianAlleleAge <- tapply(x$allele_age, x$chromEnd, quantile, 0.5)
lowerBound <- tapply(x$allele_age, x$chromEnd, quantile, 0.05)
upperBound <- tapply(x$allele_age, x$chromEnd, quantile, 0.95)

## duffy FY*O snp is at position 159174683 - commands below give estimates of age in ky
pos <- "159174683"
medianAlleleAge[[pos]]*29/1000
lowerBound[[pos]]*29/1000
upperBound[[pos]]*29/1000

