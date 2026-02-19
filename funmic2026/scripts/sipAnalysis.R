## let's look at our SIP data now that it is cleaned up and in a phyloseq object

## make sure conda is in the base environment:

conda deactivate

## go back to our metabarcoding working directory:


cd /vol/funmic/metabarcoding

## start R and get libraries:

R 


## where are you? 

getwd()

## is it where you want to be?

library(phyloseq)
library(ggplot2)

## let's get some extra functions from github repo:

download.file("https://raw.githubusercontent.com/danchurch/FunctionalMicrobiomePractical/refs/heads/main/funmic2026/scripts/sipFunctions.R",
              "sipFunctions.R", method="wget")

source("sipFunctions.R") 

## load up our phyloseq object from our metabarcoding pipeline:

load("sipPS.rda")

## how do our sequencing depths look?

plot_bar(ps)

## how might this affect our analysis?

## one very simple adjustment would be to change to proportions:
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

## a quick look at the relative abundance of the most common ASVs:

familiesPlot <- plotFamilies(ps,30)

familiesPlot 

## this takes forever. To download a copy:
ggsave("familiesBySubIso.pdf",
        device = pdf,
        plot = familiesPlot,
        width = 10,
        height = 6)


## we can try to an overview with an NMS graphic of all of our samples:
dev.new() 
NMS_braycurtis(ps.prop)

## is the stress okay?

## but that doesn't give us much information.

## let's break up our data by substrate
ps.m <- subset_samples(ps, Substrate == "M")
ps.g <- subset_samples(ps, Substrate == "G")
ps.a <- subset_samples(ps, Substrate == "A")
ps.y <- subset_samples(ps, Substrate == "Y")


## and let's look at these with a PCA approach.
## use another custom function:

dev.new()
plotPCAWithSpecies(ps.m, ptitle="Methanol PCA")

## how would you look at the other substrates: acetate, glucose and glycerol?

dev.new()
plotPCAWithSpecies(ps.a, ptitle="Acetate PCA")

dev.new()
plotPCAWithSpecies(ps.g, ptitle="Glucose PCA")

dev.new()
plotPCAWithSpecies(ps.y, ptitle="Glycerol PCA")

## the PCA and species scoring algorithms are such that the species are positioned close to 
## the sites that they are most important in distinguishing. 

## let's check one out:
## here is a function that take the following arguments:
## ASV name, the phyloseq object, and the title you want to give it:

## ASV1 looks important:
## we can graph by Fraction or Buoyant Density :

dev.new()
getFractionAbundances("ASV1",ps.m.prop,ptitle="ASV1, Methanol",fracOrBD="Fraction")

dev.new()
getFractionAbundances("ASV1",ps.m.prop,ptitle="ASV1, Methanol",fracOrBD="BD")

## an example of an "unshifted" organism: 

dev.new()

## example of two gc-rich species, no strong evidence of uptake of heavy isotope: 
getFractionAbundances(ASVname="ASV3", whichPS=ps.m.prop, ptitle="Methanol, ASV3", fracOrBD="BD")
tax_table(ps)["ASV3",]

getFractionAbundances(ASVname="ASV7", whichPS=ps.m.prop, ptitle="Methanol, ASV7", fracOrBD="BD")
tax_table(ps)["ASV7",]

## etc etc 

## remember what to do if you want to save a plot?:

methOH_ASV33 <- getFractionAbundances(ASVname="ASV33", whichPS=ps.m.prop, ptitle="glucose, ASV33", fracOrBD="BD")

ggsave("ASV33_MethOH_buoyantDensAbundances.pdf",
        device = pdf,
        plot = methOH_ASV33)

ggsave("ASV33_MethOH_buoyantDensAbundances.png",
        device = png,
        plot = methOH_ASV33)


