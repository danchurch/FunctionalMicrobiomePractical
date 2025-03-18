## let's look at our SIP data now that it is cleaned up and in a phyloseq object

## go back to our metabarcoding working directory:
cd /vol/funmic/metabarcoding

## start R and get libraries:

R 

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


#sample_data(ps.prop)$Isotope = as.factor(sample_data(ps.prop)$Isotope)


NMS_braycurtis(ps.prop)

## let's break up our data by substrate:
samp.a <- rownames(sample_data(ps)[sample_data(ps)$Substrate == "A",])
samp.g <- rownames(sample_data(ps)[sample_data(ps)$Substrate == "G",])
samp.m <- rownames(sample_data(ps)[sample_data(ps)$Substrate == "M",])
ps.a <- prune_samples(samp.a, ps)
ps.g <- prune_samples(samp.g, ps)
ps.m <- prune_samples(samp.m, ps)

## and let's look at these with a PCA approach:

dev.new()

plotPCAWithSpecies(ps=ps.m, ptitle="Methanol PCA")

dev.new()
plotPCAWithSpecies(ps=ps.a, ptitle="Acetate PCA")

dev.new()
plotPCAWithSpecies(ps=ps.g, ptitle="Glucose PCA")

## the PCA and species scoring algorithms are such that the species are positioned close to 
## the sites that they are most important in distinguishing. 

## let's check one out:
## here is a function that take the ASV name, the phyloseq object, and the title you want to give it:

dev.new()
getFractionAbundances(ASVname="ASV1", whichPS=ps.m, ptitle="methanol, ASV21")

dev.new()
getFractionAbundances(ASVname="ASV33", whichPS=ps.g, ptitle="glucose, ASV40")

## that function reports taxonomy. But we can also get it this way, with phyloseq:

tax_table(ps)["ASV21",]

## if you want to save a plot, try the pdf() or png() functions:

png(file="fractionAbundances_ASV33_glucose.png")
getFractionAbundances(ASVname="ASV33", whichPS=ps.g, ptitle="glucose, ASV33")
dev.off()



############## I have to do this to get my file, but you probably don't, just use MobaXterm!! ###############
getFile=/vol/funmic/metabarcoding/fractionAbundances_ASV33_glucose.png
putDir=/home/daniel/Documents/teaching/funmic/scratchpad/
scp -i /home/daniel/.ssh -P 30423 -r ubuntu@129.70.51.6:$getFile $putDir
#############################################################################################################

