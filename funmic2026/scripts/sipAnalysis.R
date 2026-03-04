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

download.file("https://raw.githubusercontent.com/danchurch/FunctionalMicrobiomePractical/refs/heads/main/funmic2026/scripts/sipFunctions.R", "sipFunctions.R", method="wget")

source("sipFunctions.R") 

## load up our phyloseq object from our metabarcoding pipeline:

load("sipPS.rda")

## how do our sequencing depths look?

plot_bar(ps)




## how might this affect our analysis?

## one very simple adjustment would be to change to proportions:
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

## a quick look at the relative abundance of the most common ASVs:

familiesPlot <- plotFamilies(ps.prop,30)

#familiesPlot ## but wait!! takes a really long time.  

## since this takes forever. To get a copy to download:
ggsave("familiesBySubIso.pdf",
        device = pdf,
        plot = familiesPlot,
        width = 10,
        height = 6)



## we can try to an overview with an NMS graphic of all of our samples:
#dev.new() 

nmsAllsamples <- NMS_braycurtis(ps.prop)

nmsAllsamples

#ggsave("nmsAllsamples.pdf",
#        device = pdf,
#        plot = nmsAllsamples,
#        width = 10,
#        height = 6)

## is the stress okay?

## but that doesn't give us much information.

## let's break up our data by substrate
ps.m <- subset_samples(ps, Substrate == "M")
ps.g <- subset_samples(ps, Substrate == "G")

ps.a <- subset_samples(ps, Substrate == "A")

ps.y <- subset_samples(ps, Substrate == "Y")


## and let's look at these with a PCA approach.
## use another custom function: "plotPCAWithSpecies"



dev.new()

acetatePCA <- plotPCAWithSpecies(ps.a, ptitle="Acetate PCA")

acetatePCA 

ggsave("acetatePCA.pdf",
        device = pdf,
        plot = acetatePCA,
        width = 10,
        height = 6)



## how would you look at the other substrates: acetate, glucose and glycerol?


methanolPCA <- plotPCAWithSpecies(ps.m, ptitle="Methanol PCA")
glucosePCA <- plotPCAWithSpecies(ps.g, ptitle="Glucose PCA")
glycerolPCA <- plotPCAWithSpecies(ps.y, ptitle="Glycerol PCA")

#ggsave("glycerolPCA.pdf",
#        device = pdf,
#        plot = glycerolPCA,
#        width = 10,
#        height = 6)

## the PCA and species scoring algorithms are such that the species are positioned close to 
## the sites that they are most important in distinguishing. 




## remember, we haven't transformed our abundances of our phyloseq objects by substrate
## how would we do this?

ps.m.prop <- transform_sample_counts(ps.m, function(otu) otu/sum(otu))
ps.g.prop <- transform_sample_counts(ps.g, function(otu) otu/sum(otu))

ps.a.prop <- transform_sample_counts(ps.a, function(otu) otu/sum(otu))

ps.y.prop <- transform_sample_counts(ps.y, function(otu) otu/sum(otu))


## let's check one out:
## here is a function that take the following arguments:
## ASV name, the phyloseq object, and the title you want to give it:

## ASV1 looks important:
## we can graph by Fraction or Buoyant Density :


ASV1Methanol_fraction <- getFractionAbundances("ASV1",ps.m.prop,ptitle="ASV1, Methanol",fracOrBD="Fraction")

dev.new()
ASV1Methanol_fraction


ASV1Methanol_BD  <- getFractionAbundances("ASV1",ps.m.prop,ptitle="ASV1, Methanol",fracOrBD="BD")

dev.new()
ASV1Methanol_BD

ASV100Acetate_BD  <- getFractionAbundances("ASV100",ps.a.prop,ptitle="ASV100, Acetate",fracOrBD="BD")
ASV455Acetate_BD <- getFractionAbundances("ASV455",ps.a.prop,ptitle="ASV455, Acetate",fracOrBD="BD")
ASV31Acetate_BD  <- getFractionAbundances("ASV31",ps.a.prop,ptitle="ASV31, Acetate",fracOrBD="BD")
ggsave("ASV31Acetate_BD.pdf",
        device = pdf,
        plot = ASV31Acetate_BD,
        width = 10,
        height = 6)

#getFile="/vol/funmic/metabarcoding/acetatePCA.pdf"
getFile="/vol/funmic/metabarcoding/ASV*Acetate_BD.pdf"
putDir="/home/daniel/Documents/teaching/functionalMicrobiomes/scratchpad/"
rsync -auvn \
      --progress \
      -e "ssh -p 31993" \
      ubuntu@129.70.51.6:$getFile $putDir

## this looks better. 

## an example of an "unshifted" organism: 

dev.new()

## example of two gc-rich species, no strong evidence of uptake of heavy isotope: 

getFractionAbundances(ASVname="ASV3", whichPS=ps.m.prop, ptitle="Methanol, ASV3", fracOrBD="BD")

tax_table(ps)["ASV3",]

getFractionAbundances(ASVname="ASV7", whichPS=ps.m.prop, ptitle="Methanol, ASV7", fracOrBD="BD")
tax_table(ps)["ASV7",]

## etc etc 


## getting taxonomy:

tax_table(ps)["ASV1"]



### getting out spreadsheets of our data for tomorrow's "SIP amplicon data mining" session:


