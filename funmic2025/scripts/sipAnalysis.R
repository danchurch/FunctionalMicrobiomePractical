## let's look at our SIP data now that it is cleaned up and in a phyloseq object

## go back to our metabarcoding working directory:
cd /vol/funmic/metabarcoding

## start R and get libraries:

R 

library(phyloseq)
source("/vol/funmic/scripts/sipFunctions.R")

## load up our phyloseq object from our metabarcoding pipeline:

load(
