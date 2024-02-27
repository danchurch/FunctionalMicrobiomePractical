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
			--out-prefix ${SAMPLE} > /vol/funMicStorage/Kelp_logs/${SAMPLE}_megahitlog.txt
			
	done
	
## Time with 7 cores: 53 min for k21, 30 min for subsequent iterations, 9h 40 min per metagenome
## Simplify contig names down to simple IDs, removing everything that follows the space in the header
sed -i 's/\ .*//g' assemblies/kelp/*/*.contigs.fa
for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
	do
		sed -i 's/>/>'$ASSEMBLY'_/g' assemblies/kelp/${ASSEMBLY}/${ASSEMBLY}.contigs.fa
	done

## okay, let's look at it with quast

conda activate quast

mkdir /vol/funMicStorage/assemblyQC
ASSEMBLY_DIRECTORY=/vol/funMicStorage/assemblies/kelp/
OUTPUTDIRECTORY=/vol/funMicStorage/assemblyQC/

for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
	do 
		quast -t 12 \
			-o ${OUTPUTDIRECTORY}/${ASSEMBLY}/ \
			${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
			> /vol/funMicStorage/Kelp_logs/${ASSEMBLY}_quastLog.txt
	done

conda deactivate	
## trivial amount of time and memory used.

####### Read mapping for differential coverage #########
conda activate RawReadProcessing

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
				ref=${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
				in=${READ_DIRECTORY}/${SAMPLE}_1.fastq.gz \
				in2=${READ_DIRECTORY}/${SAMPLE}_2.fastq.gz \
				outm=${OUTPUTDIRECTORY}/${SAMPLE}_to_${ASSEMBLY}.sam \
				bamscript=${OUTPUTDIRECTORY}/${SAMPLE}_to_${ASSEMBLY}.sh \
				> /vol/funMicStorage/Kelp_logs/${SAMPLE}_to_${ASSEMBLY}_bbmaplog.txt
				
				bash ${OUTPUTDIRECTORY}/${SAMPLE}_to_${ASSEMBLY}.sh
				rm ${OUTPUTDIRECTORY}/*.sam
			done
	done

conda deactivate
##Time: 3h 20 min

####### binning ######

## let's try metabat on a couple different settings

conda activate binning

mkdir /vol/funMicStorage/binning

#defining some variables. DATA_FOLDER should be an absolute path to where the data for the project is stored, with fastq reads being in a subfolder called "reads".
#PROJECTNAME is a prefix which was used for the assembly etc.
#METAWATTPATH should be the absolute path to the folder with the Metawatt jar file. METAWATT_DATABASES is the absolute path to the Metawatt databases
#CONCOCTPATH should be the absolute path to the concoct folder containing subfolders "bin" and "scripts"

MAPPING_DIRECTORY=/vol/funMicStorage/mapping
ASSEMBLY_DIRECTORY=/vol/funMicStorage/assemblies/kelp
BINNING_DIRECTORY=/vol/funMicStorage/binning


#Binning of the metagenomes co-assembly 

#This script is part of the MetaBAT module and reads out the bam files to generate a coverage table (incl. coverage variation) for the contigs

for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
	do  
		mkdir ${BINNING_DIRECTORY}/${ASSEMBLY}
		mkdir ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT
		jgi_summarize_bam_contig_depths \
			--outputDepth ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_coverage_Depths.txt \
			--referenceFasta  ${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
			${MAPPING_DIRECTORY}/*_to_${ASSEMBLY}_sorted.bam
    done
    
#run MetaBAT
for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
	do     
	metabat2 \
        -i ${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
        -a ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_coverage_Depths.txt \
        -o ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_MetaBAT \
        -t 12
	done

## Time neglidgible (few minutes)

## perform binning with vamb
for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
	do     
	vamb \
        --fasta ${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
        --jgi ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_coverage_Depths.txt \
        --outdir ${BINNING_DIRECTORY}/${ASSEMBLY}/vamb/ \
		--minfasta 200000
	done

## Time: 2 h 7 min.

## perform binning with CONCOCT
conda deactivate
conda activate concoct

MAPPING_DIRECTORY=/vol/funMicStorage/mapping
ASSEMBLY_DIRECTORY=/vol/funMicStorage/assemblies/kelp
BINNING_DIRECTORY=/vol/funMicStorage/binning

## create contig coverage tables for CONCOCT from the MetaBAT coverage tables.
## The first column with coverage is actually column 4.
## Start from column 4, pick out every second column (CONCOCT doesn't need the coverage variation columns, only the coverage itself) until you got all the samples. 
## Also remove the first line with the headers.

for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
	do 
		mkdir ${BINNING_DIRECTORY}/${ASSEMBLY}/concoct
		cut -f 1,4,6,8 ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_coverage_Depths.txt > ${BINNING_DIRECTORY}/${ASSEMBLY}/concoct/${ASSEMBLY}_coverage_Depths.txt
	done

for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
	do     
	concoct \
		-t 12 \
        --composition_file ${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
        --coverage_file ${BINNING_DIRECTORY}/${ASSEMBLY}/concoct/${ASSEMBLY}_coverage_Depths.txt \
		-r 250 \
        -b ${BINNING_DIRECTORY}/${ASSEMBLY}/concoct/${ASSEMBLY}_CONCOCT \
        --no_original_data
	done

## Time: 4 h

## Choosing best bins from the results of different binning tools for each assembly with DAStool. DAStool judges bin size, completeness and contamination values
## generate DAS_Tool input files: Tables with two columns assigning each contig to a putative genomic bin. For each assembly we need such tables from every binning tool.
## Attribution of contigs to bins according to MetaBAT
for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
  do 
	mkdir ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool
	  for BIN in $(ls ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_MetaBAT.*.fa | xargs -n 1 basename | sed 's/\.fa//g'); 
		do 
			grep '^>' ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/"$BIN".fa | sed "s/$/\t$BIN/g" | sed 's/>//g' \
					>> ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool/${ASSEMBLY}_MetaBAT_binning.txt 
		done
  done

## Attribution of contigs to bins according to vamb
for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
  do 
	  for BIN in $(ls ${BINNING_DIRECTORY}/${ASSEMBLY}/vamb/bins/*.fna | xargs -n 1 basename | sed 's/\.fna//g'); 
		do 
			grep '^>' ${BINNING_DIRECTORY}/${ASSEMBLY}/vamb/bins/${BIN}.fna | sed "s/$/\t$ASSEMBLY\_vamb$BIN/g" | sed 's/>//g' \
					>> ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool/${ASSEMBLY}_vamb_binning.txt 
		done
  done


## Attribution of contigs to bins according to CONCOCT
for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
  do 
	sed '1d' ${BINNING_DIRECTORY}/${ASSEMBLY}/concoct/${ASSEMBLY}_CONCOCT_clustering_gt1000.csv | \
	perl -0pe "s/\,/\t$ASSEMBLY\_CONCOCT/g" \
	> ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool/${ASSEMBLY}_CONCOCT_binning.txt
  done

conda deactivate

## Run DAS_Tool to generate consensus bins from the results of different binners
conda activate binning

ASSEMBLY_DIRECTORY=/vol/funMicStorage/assemblies/kelp
BINNING_DIRECTORY=/vol/funMicStorage/binning

for ASSEMBLY in $(ls /vol/funMicStorage/assemblies/kelp/ | xargs -n 1 basename);
  do 
	BINNING_FILES=$(ls ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool/*_binning.txt | \
					perl -0pe 's/\n/\,/g' | \
					sed 's/\,$//g')

	BINNERS=$(ls ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool/*_binning.txt | xargs -n 1 basename | \
			sed "s/${ASSEMBLY}\_//g" | \
			sed 's/\_binning.txt//g' | \
			perl -0pe 's/\n/\,/g' | \
			sed 's/\,$//g')

	DAS_Tool \
		-i "$BINNING_FILES" \
		-l "$BINNERS" \
		-c ${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
		-o ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool/${ASSEMBLY} \
		--score_threshold 0.25 \
		-t 12 \
		--write_bins
  done

## Time 7 min

conda deactivate


## Dereplicate bins from different samples
## As the community is very similar between the samples, same genomes will be obtained from each of the assemblies
## dRep identifies "identical" bins and keeps the bins with the best metrics (size, completeness, contamination)
## In the end we should obtain one optimized set of non-redundant bins from the three different samples.

conda activate drep

mkdir /vol/funMicStorage/binning/drep
mkdir /vol/funMicStorage/binning/drep/all_bins
mkdir /vol/funMicStorage/binning/drep/dereplicated_bins

BINNING_DIRECTORY=/vol/funMicStorage/binning
ALL_BINS_DIRECTORY=/vol/funMicStorage/binning/drep/all_bins
DEREP_BINS_DIRECTORY=/vol/funMicStorage/binning/drep/dereplicated_bins

scp ${BINNING_DIRECTORY}/ERR3801*/DAStool/*_bins/*.fa ${ALL_BINS_DIRECTORY}/

export OMP_NUM_THREADS=12
export NUMEXPR_MAX_THREADS=12

dRep dereplicate \
    ${DEREP_BINS_DIRECTORY} \
    -p 12 \
    -g ${ALL_BINS_DIRECTORY}/*.fa \
    -comp 40 \
    -l 500000 \
    -comW 1 \
    -sizeW 1 \
    -conW 5 \
    -strW 1 \
    -N50W 0.5 \
    --run_tertiary_clustering

## Time 1 h 20 min

conda deactivate

## Check MAG completeness and contamination with the new checkm2. Results from CheckM 1.2.2 should be available from dRep

conda activate magQC

BINNING_DIRECTORY=/vol/funMicStorage/binning
ALL_BINS_DIRECTORY=/vol/funMicStorage/binning/drep/all_bins

DEREP_BINS_DIRECTORY=/vol/funMicStorage/binning/drep/dereplicated_bins/dereplicated_genomes

checkm2 predict \
	--allmodels \
	--threads 12 \
	--remove_intermediates \
	--input ${DEREP_BINS_DIRECTORY} \
	--extension .fa \
	--output-directory ${BINNING_DIRECTORY}/dereplicated_bins_CheckM
	
## Time 6 min
conda deactivate

##  Classify the bins taxonomically with GTDB-Tk
conda activate gtdbtk

BINNING_DIRECTORY=/vol/funMicStorage/binning
DEREP_BINS_DIRECTORY=/vol/funMicStorage/binning/drep/dereplicated_bins/dereplicated_genomes


gtdbtk classify_wf \
	--cpus 12 \
	--pplacer_cpus 12 \
	--genome_dir ${DEREP_BINS_DIRECTORY} \
	-x fa \
	--out_dir ${BINNING_DIRECTORY}/dereplicated_bins_GTDB
	
## Time 4h 45 min


####-----------------------------------------------------------------------------------------------------------------------
#### Generating an Anvi'O database (right now Anvi'O  interactive mode is not working)
####-----------------------------------------------------------------------------------------------------------------------

mkdir /vol/funMicStorage/Anvio_kelp
mkdir /vol/funMicStorage/Anvio_kelp/input

ANVIO_DIRECTORY=/vol/funMicStorage/Anvio_kelp
DEREP_BINS_DIRECTORY=/vol/funMicStorage/binning/drep/dereplicated_bins/dereplicated_genomes
READ_DIRECTORY=/vol/funMicStorage/datasets/kelp

## Combine all bins into one FASTA FILE

cat ${DEREP_BINS_DIRECTORY}/*.fa > ${ANVIO_DIRECTORY}/input/All_bins.fasta

## Map the reads from all samples to the combined fasta file

conda activate RawReadProcessing

ANVIO_DIRECTORY=/vol/funMicStorage/Anvio_kelp
DEREP_BINS_DIRECTORY=/vol/funMicStorage/binning/drep/dereplicated_bins/dereplicated_genomes
READ_DIRECTORY=/vol/funMicStorage/datasets/kelp

for SAMPLE in $(ls ${READ_DIRECTORY}/*_1.fastq.gz | xargs -n 1 basename | sed 's/_1.fastq.gz//g');
	do 
		bbmap.sh \
		threads=12 \
		minid=97 \
		idfilter=95 \
		ref=${ANVIO_DIRECTORY}/input/All_bins.fasta \
		in=${READ_DIRECTORY}/${SAMPLE}_1.fastq.gz \
		in2=${READ_DIRECTORY}/${SAMPLE}_2.fastq.gz \
		outm=${ANVIO_DIRECTORY}/input/${SAMPLE}.sam \
		bamscript=${ANVIO_DIRECTORY}/input/${SAMPLE}.sh
				
		bash ${ANVIO_DIRECTORY}/input/${SAMPLE}.sh
		rm ${ANVIO_DIRECTORY}/*.sam
	done

conda deactivate
conda activate anvio-8

ANVIO_DIRECTORY=/vol/funMicStorage/Anvio_kelp

## Create the contig database
anvi-gen-contigs-database \
	-T 12 \
	-L 0 \
	-n Kelp \
	-f ${ANVIO_DIRECTORY}/input/All_bins.fasta \
	-o ${ANVIO_DIRECTORY}/Kelp.db

## Create coverage ("profile") dadabases for each sample and merge them into one profile database

for SAMPLE in $(ls ${ANVIO_DIRECTORY}/input/*_sorted.bam | xargs -n 1 basename | sed 's/_sorted.bam//g'); 
	do 
		anvi-profile \
		-T 12 \
		-M 1000 \
		-c ${ANVIO_DIRECTORY}/Kelp.db \
		-o ${ANVIO_DIRECTORY}/Individual_Profiles/${SAMPLE} \
		--sample-name ${SAMPLE} \
		-i ${ANVIO_DIRECTORY}/input/${SAMPLE}_sorted.bam 
	done
	
anvi-merge \
	-c ${ANVIO_DIRECTORY}/Kelp.db \
	-o ${ANVIO_DIRECTORY}/PROFILE \
	${ANVIO_DIRECTORY}/Individual_Profiles/*/PROFILE.db

anvi-run-hmms \
	-T 12 \
	-c ${ANVIO_DIRECTORY}/Kelp.db	
	
anvi-run-scg-taxonomy \
	-c ${ANVIO_DIRECTORY}/Kelp.db \
	-P 12

DEREP_BINS_DIRECTORY=/vol/funMicStorage/binning/drep/dereplicated_bins/dereplicated_genomes

for BIN in $(ls ${DEREP_BINS_DIRECTORY}/*.fa | xargs -n 1 basename | sed 's/\.fa//g'); 
		do 
			grep '^>' ${DEREP_BINS_DIRECTORY}/"$BIN".fa | sed "s/$/\t$BIN/g" | sed 's/>//g' \
					>> ${ANVIO_DIRECTORY}/input/Dereplicated_bins_collection.txt 
		done

## Export the translated amino acid sequences for gene calls, to run a taxonomic classification

anvi-get-sequences-for-gene-calls \
	--get-aa-sequences \
	-c ${ANVIO_DIRECTORY}/Kelp.db \
	-o ${ANVIO_DIRECTORY}/Kelp.faa



		
## Anvi'O doesn't accept . in the MAG names, which is the case for MetaBAT MAGs (MetaBAT.1, MetaBAT.12 etc.)
## Search and replace dots with nothing (aka remove the dots).

sed -i 's/\.//g' ${ANVIO_DIRECTORY}/input/Dereplicated_bins_collection.txt

## Import the MAG collection twice under a different name: one will contain the original MAGs, the other will be for MAG refinement
anvi-import-collection \
	-c ${ANVIO_DIRECTORY}/Kelp.db \
	-p ${ANVIO_DIRECTORY}/PROFILE/PROFILE.db \
	--contigs-mode \
	-C Dereplicated_bins \
	${ANVIO_DIRECTORY}/input/Dereplicated_bins_collection.txt

anvi-import-collection \
	-c ${ANVIO_DIRECTORY}/Kelp.db \
	-p ${ANVIO_DIRECTORY}/PROFILE/PROFILE.db \
	--contigs-mode \
	-C Refined_bins \
	${ANVIO_DIRECTORY}/input/Dereplicated_bins_collection.txt

## Find a MAG with considerable contamination and try to refine it

