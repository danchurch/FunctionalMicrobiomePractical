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

## The Illumina reads from the white kelp biofims from Helgoland is here: /vol/funmic/datasets/kelpBiofilm
## We store this location under the variable READ_DIRECTORY
READ_DIRECTORY=/vol/funmic/datasets/kelpBiofilm

#Checking read quality with FastQC for all reads. FastQC can also be run with an interactive GUI by typing 'fastqc' into the command line
mkdir /vol/funmic/Kelp/rawReadQC
OUTPUTDIRECTORY=/vol/funmic/Kelp/rawReadQC

## The variable FILE is read from listing all files ending with fastq.gz in the folder containing the data.
## xargs ensures only the file names and not the entire path is listed
## sed performs a search and replace operation replacing .fastq.gz with nothing, so that only the sample names are read into the vaiable and the file type extension is not dragged along

conda activate qualityControlRawSequences

for FILE in $(ls ${READ_DIRECTORY}/*.fastq.gz | xargs -n 1 basename | sed 's/.fastq.gz//g');
	do
		fastqc \
		-t 12 \
		-o ${OUTPUTDIRECTORY} \
		${READ_DIRECTORY}/${FILE}.fastq.gz > /vol/funmic/Kelp/logs/${FILE}.fastqcLog.txt
	done

conda deactivate

## The quality statistics are summarized in the fastqc.html files for each read file
## Can be opened in a browser, for example by navigating to the containing folder in MobaXterm file browser (left side), right-clicking the html and select "open with " e.g. Chrome


####### assembly #####

conda activate assembly

mkdir -p /vol/funmic/Kelp/assemblies

READ_DIRECTORY=/vol/funmic/datasets/kelpBiofilm
OUTPUTDIRECTORY=/vol/funmic/Kelp/assemblies

## Iterative assembly with different k-mer sizes. K-mer size must be odd. K-mer sizes range from 21 and up to the full read length in increments of 10 bp
## At the end contigs shorter than 1000 bp are discarded from the final assembly, as they cannot be reliably attributed to genomes or annotated in a meaningful way

for SAMPLE in $(ls ${READ_DIRECTORY}/*_1.fastq.gz | xargs -n 1 basename | sed 's/_1.fastq.gz//g');
	do
		megahit \
			-1 ${READ_DIRECTORY}/${SAMPLE}_1.fastq.gz \
			-2 ${READ_DIRECTORY}/${SAMPLE}_2.fastq.gz \
			-o ${OUTPUTDIRECTORY}/${SAMPLE} \
			--k-min 21 \
			--k-max 249 \
			--k-step 10 \
			-m 1 \
			-t 12 \
			--min-contig-len 1000 \
			--out-prefix ${SAMPLE} > /vol/funmic/Kelp/logs/${SAMPLE}_megahitlog.txt

	done



