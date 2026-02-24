#!/bin/bash

#### Metagenome analysis #######


## In this course we are analyzing DNA shotgun sequencing data from the following publication:
## paper: Van Erk et al. 2021, doi.org/mgt4
## dataset accession number on NCBI/ENA: PRJEB36085

ssh deNBI

## First, create a folder  for all the new files that you will generate during the analysis.
## We will keep it organized and structured to easily find the output produced by the different tools

# Make a new directory "Kelp" under /vol/funmic/
mkdir /vol/funmic/Kelp

# Make a sub directory to store log files of some programs
mkdir /vol/funmic/Kelp/logs

