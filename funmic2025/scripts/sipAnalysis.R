## let's look at our SIP data now that it is cleaned up and in a phyloseq object

## go back to our metabarcoding working directory:

cd /vol/funmic/metabarcoding

## start R and get libraries:

R 

## just make sure we are in the right place
setwd("/vol/funmic/metabarcoding") 

library(phyloseq)
library(ggplot2)
source("/vol/funmic/scripts/sipFunctions.R") ## some extra functions

## load up our phyloseq object from our metabarcoding pipeline:

load("sipPS.rda")

## how do our sequencing depths look?

plot_bar(ps)

## how might this affect our analysis?

## one very simple adjustment would be to change to proportions:
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

## a quick look at the relative abundance of the most common ASVs:

plotFamilies(30)

dev.new() ## start a new plotter, so we don't clobber the old figure

NMS_braycurtis(ps.prop)


## let's break up our data by substrate, making abundances in each proportional:

## methanol
samp.m <- rownames(sample_data(ps)[sample_data(ps)$Substrate == "M",])
ps.m <- prune_samples(samp.m, ps)
ps.m.prop <- transform_sample_counts(ps.m, function(otu) otu/sum(otu))

## glucose
samp.g <- rownames(sample_data(ps)[sample_data(ps)$Substrate == "G",])
ps.g <- prune_samples(samp.g, ps)
ps.g.prop <- transform_sample_counts(ps.g, function(otu) otu/sum(otu))

## acetate
samp.a <- rownames(sample_data(ps)[sample_data(ps)$Substrate == "A",])
ps.a <- prune_samples(samp.a, ps)
ps.a.prop <- transform_sample_counts(ps.a, function(otu) otu/sum(otu))


## and let's look at these with a PCA approach:

dev.new()

plotPCAWithSpecies(ps.m.prop, ptitle="Methanol PCA")

dev.new()
plotPCAWithSpecies(ps.a.prop, ptitle="Acetate PCA")

dev.new()
plotPCAWithSpecies(ps.g.prop, ptitle="glucose PCA")

## the PCA and species scoring algorithms are such that the species are positioned close to 
## the sites that they are most important in distinguishing. 

## let's check one out:
## here is a function that take the ASV name, the phyloseq object, and the title you want to give it:


## plots from yesterday looked like this:

dev.new()

getFractionAbundances("ASV1",ps.m.prop,ptitle="ASV1, Methanol",fracOrBD="Fraction")

## we added the ability to graph by Buoyant Density:

dev.new()
getFractionAbundances("ASV1",ps.m.prop,ptitle="ASV1, Methanol",fracOrBD="BD")

dev.new()

getFractionAbundances(ASVname="ASV19", whichPS=ps.g.prop, ptitle="glucose, ASV19", fracOrBD="BD")

getFractionAbundances("ASV1",ps.m.prop,ptitle="ASV1, Methanol",fracOrBD="BD")


## that function reports taxonomy. But we can also get it this way, with phyloseq:

tax_table(ps)["ASV8",]

## remember that if you want to save a plot, try the pdf() or png() functions:

png(file="fractionAbundances_ASV33_glucose.png")
getFractionAbundances(ASVname="ASV33", whichPS=ps.g.prop, ptitle="glucose, ASV33", fracOrBD="BD")
dev.off()


############## I have to do this to get my file, but you probably don't, just use MobaXterm!! ###############
getFile=/vol/funmic/metabarcoding/fractionAbundances_ASV33_glucose.png
putDir=/home/daniel/Documents/teaching/funmic/scratchpad/
scp -i /home/daniel/.ssh -P 30423 -r ubuntu@129.70.51.6:$getFile $putDir
#############################################################################################################

