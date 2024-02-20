#!/bin/bash

#### Metagenome analysis #######

## I. Checking the Illumina sequence reads: quality of the data? Microbial composituion of the raw data: What did get sequenced?

conda activate RawReadProcessing

## The Illumina reads from the white kelp biofims from Helgoland is here:
READ_DIRECTORY=/vol/funMicStorage/datasets/kelp

#paper: Van Erk et al. 2021, doi.org/mgt4
#dataset accession number on NCBI/ENA: PRJEB36085

#Checking read quality with FastQC for all reads. FastQC can also be run with an interactive GUI by typing 'fastqc' into the command line
mkdir /vol/funMicStorage/rawReadQC
OUTPUTDIRECTORY=/vol/funMicStorage/rawReadQC

## The variable FILE is read from listing all files ending with fastq.gz in the folder containing the data. 
## xargs ensures only the file names and not the entire path is listed
## sed performs a search and replace operation replacing .fastq.gz with nothing, so that only the sample names are read into the vaiable and the file type extension is not dragged along

for FILE in $(ls ${READ_DIRECTORY}/*.fastq.gz | xargs -n 1 basename | sed 's/.fastq.gz//g');
	do 
		fastqc \
		-t 7 \
		-o ${OUTPUTDIRECTORY} \
		${READ_DIRECTORY}/${FILE}.fastq.gz > /vol/funMicStorage/Kelp_logs/${FILE}.fastqcLog.txt
	done

#Time: 22 min

## The quality statistics are summarized in the fastqc.html files for each read file
## Can be opened in a browser, for example by navigating to the containing folder in MobaXterm file browser (left side), right-clicking the html and select "open with " e.g. Chrome 

##### Analyzing taxonomic composition of the datasets with phyloflash: mapping of reads to small subunit ribosomal RNA database 
## and targeted assembly SSU rRNA sequences #####

phyloFlash.pl -man

## Defining the location of the SSU rRNA database 
##(SILVA non-redundnat reference sequences, edition 138.1, pre-processed and formatted for phyloFlash)
phyloDB=/vol/funMicStorage/databases/phyloFlash/138.1
READ_DIRECTORY=/vol/funMicStorage/datasets/kelp
OUTPUTDIRECTORY=/vol/funMicStorage/phyloFlashOut/

## The variable SAMPLE is read from listing all files ending with _1.fastq.gz in the folder containing the data.
## As reads for each sample are split into two files, we should only list one of the pair in order not to have samples listed twice.
 
for SAMPLE in $(ls ${READ_DIRECTORY}/*_1.fastq.gz | xargs -n 1 basename | sed 's/_1.fastq.gz//g');
	do 
		mkdir ${OUTPUTDIRECTORY}/${SAMPLE}
		phyloFlash.pl \
			-lib PhyloFlash_Kelp_${SAMPLE} \
			-read1 ${READ_DIRECTORY}/${SAMPLE}_1.fastq.gz \
			-read2 ${READ_DIRECTORY}/${SAMPLE}_2.fastq.gz \
			-readlength 250 \
			-clusterid 98 \
			-taxlevel 7 \
			-dbhome $phyloDB \
			-CPUs 7 > /vol/funMicStorage/Kelp_logs/${SAMPLE}_phyloflashlog.txt
		
		mv PhyloFlash_Kelp_${SAMPLE}* ${OUTPUTDIRECTORY}/${SAMPLE}
	done
	
##Time: 45 min ##
conda deactivate

####### assembly #####

conda activate assembly

mkdir -p /vol/funMicStorage/assemblies/kelp/

READ_DIRECTORY=/vol/funMicStorage/datasets/kelp
OUTPUTDIRECTORY=/vol/funMicStorage/assemblies/kelp/

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
			-t 7 \
			--min-contig-len 1000 \
			--out-prefix ${SAMPLE} > /vol/funMicStorage/Kelp_logs/${SAMPLE}_megahitlog.txt
			
	done
	
##Time: 11:24:51 - (53 min for k21, 30 min for subsequent iterations, 9h 40 min per metagenome)

## okay, let's look at it with quast

conda activate quast

mkdir /vol/funMicStorage/assemblyQC
ASSEMBLY_DIRECTORY=/vol/funMicStorage/assemblies/kelp/
OUTPUTDIRECTORY=/vol/funMicStorage/assemblyQC/

for SAMPLE in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
	do 
		quast -t 7 \
			-o ${OUTPUTDIRECTORY}/${SAMPLE}/ \
			${ASSEMBLY_DIRECTORY}/${SAMPLE}/${SAMPLE}.contigs.fa \
			> /vol/funMicStorage/Kelp_logs/${SAMPLE}_quastLog.txt
	done
	
## trivial amount of time and memory used.

####### Read mapping for differential coverage #########
mkdir /vol/funMicStorage/mapping
ASSEMBLY_DIRECTORY=/vol/funMicStorage/assemblies/kelp/
READ_DIRECTORY=/vol/funMicStorage/datasets/kelp
OUTPUTDIRECTORY=/vol/funMicStorage/mapping

## each of the three read sets is mapped to each of the three assemblies: 9 mapping operations in total
## minid specifies the approximate nucleotide sequence identity to aim for when mapping reads. 
## Real identity of the alignment may be slightly higher or lower.
## idfilter parameter sets the real cut-off for reported alignments/mappings. Reads that map with an nucleotide identity less than 95% won't be reported.
## 95% is the average nucleotide identity between microbial genomes of the same species (e.g. different strains of the same species).

for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
	do
		for SAMPLE in $(ls ${READ_DIRECTORY}/*_1.fastq.gz | xargs -n 1 basename | sed 's/_1.fastq.gz//g');
			do 
				bbmap.sh \
				threads=12 \
				minid=97 \
				idfilter=95 \
				ref=${ASSEMBLY_DIRECTORY}/${SAMPLE}/${SAMPLE}.contigs.fa \
				in=${READ_DIRECTORY}/${SAMPLE}_1.fastq.gz \
				in2=${READ_DIRECTORY}/${SAMPLE}_2.fastq.gz \
				outm=${OUTPUTDIRECTORY}/${SAMPLE}_to_${ASSEMBLY}.sam \
				bamscript=${OUTPUTDIRECTORY}/${SAMPLE}_to_${ASSEMBLY}.sh \
				> /vol/funMicStorage/Kelp_logs/${SAMPLE}_to_${ASSEMBLY}_bbmaplog.txt
				
				bash ${OUTPUTDIRECTORY}/${SAMPLE}_to_${ASSEMBLY}.sh
				rm ${OUTPUTDIRECTORY}/*.sam
			done
	done


##Time: 3h 20 min

####### binning ######

## let's try metabat on a couple different settings
## and VAMB on default

## VAMB ##

conda activate binning

mkdir /vol/funMicStorage/binning

cd /vol/funMicStorage/binning

## VAMB requires a different format for the alignment of
## read than we used above for metabat

conda deactivate
conda activate vamb

vamb --outdir vambOut \
  --fasta illumcatalogue.fna.gz \
  --bamfiles illumReadsAligned2Contigs.bam \
  -t 8 \
  -o C \
  --minfasta 200000

## -t 8 is for very simple communities, this number can be much higher