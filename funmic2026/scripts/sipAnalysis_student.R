## this is a continuation of our metabarcoding analysis from last week. 

## where is this file

find . -name sipPS.rda

R



## we need some extra, custom functions for handling our SIP data:

download.file("https://raw.githubusercontent.com/danchurch/FunctionalMicrobiomePractical/refs/heads/main/funmic2026/scripts/sipFunctions.R", "sipFunctions.R", method="wget")

source("sipFunctions.R")

library(phyloseq)
library(ggplot2)
load("sipPS.rda")





## let's look at our overall read abundances:



plot_bar(ps)





## do you think these differences in read abundances might cause problems?










## one very simple adjustment would be to change to proportions:
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))


ps 

ps.prop 




## a quick look at abundances in each sample by family:

## a quick look at the relative abundance of the most common ASVs:




familiesPlot <- plotFamilies(ps.prop,30)


## plotting interactively takes forever. To get a copy to download:
ggsave("familiesBySubIso.pdf",
        device = pdf,
        plot = familiesPlot,
        width = 10,
        height = 6)







## as a first step, let's ordinate all of our samples into a single
## graphic:

nmsAllsamples <- NMS_braycurtis(ps.prop) ## what argument goes in here???


nmsAllsamples



## you can try to view it in the plotter or save the ggplot graphic as we
## did above with ggsave()





## do you "trust" this NMS? why or why not?







## let's break up our phyloseq object by substrates

## example, methanol phyloseq

ps.m <- subset_samples(ps, Substrate == "M")
ps.g <- subset_samples(ps, Substrate == "G")
ps.a <- subset_samples(ps, Substrate == "A")
ps.y <- subset_samples(ps, Substrate == "Y")



## how would you subset your phyloseq to other substrates?:









## and let's look at these with a PCA approach.
## use another custom function: "plotPCAWithSpecies"

## these PCAs don't need transformed data, it is already
## built into the functions. So just use your new 
## subsetted phyloseq objects. z.B.:

acetatePCA <- plotPCAWithSpecies(ps.a, ptitle="Acetate PCA")


## how would you plot the ordinations of other substrates?









## the PCA and species scoring algorithms are such that the species are positioned close to
## the sites that they are most important in distinguishing.


## we can graph patterns of abundance in individual species with another
## custom function, "getFractionAbundances"


## however, this function doesn't have transformations built into it.
## remember, we haven't transformed our abundances of our phyloseq objects by substrate
## how would we do this?

ps.m.prop <- transform_sample_counts(ps.m, function(otu) otu/sum(otu))

## and the others?:






## the "getFractionAbundances" function takes some arguments:
## 1. which ASV you want to look at, 
## 2. which phyloseq object (=which substrate)  
## 3. give it a title
## 4. tell it whether you want X-axis as fraction number or buoyant density


## ASV 1 looks important in the methanol-enriched samples...

ASV1Methanol_fraction <- getFractionAbundances("ASV1",ps.m.prop,ptitle="ASV1, Methanol",fracOrBD="Fraction")

ASV1Methanol_fraction


dev.new()

ASV1Methanol_fraction


## with BD instead:
ASV1Methanol_BD  <- getFractionAbundances("ASV1",ps.m.prop,ptitle="ASV1, Methanol",fracOrBD="BD")

dev.new()

ASV1Methanol_BD  




## this function returns the taxonomy of this ASV, but it can get lost. 
## so to get this taxonomy of this organism:

tax_table(ps)["ASV1"]


## what do unshifted species look like?

tax_table(ps)["ASV3",]

dev.new()

getFractionAbundances(ASVname="ASV3", whichPS=ps.m.prop, ptitle="Methanol, ASV3", fracOrBD="BD")



