#!/bin/bash

#### Metagenome analysis #######


## In this course we are analyzing DNA shotgun sequencing data from the following publication:
## paper: Van Erk et al. 2021, doi.org/mgt4
## dataset accession number on NCBI/ENA: PRJEB36085

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

#Time: 22 min

## The quality statistics are summarized in the fastqc.html files for each read file
## Can be opened in a browser, for example by navigating to the containing folder in MobaXterm file browser (left side), right-clicking the html and select "open with " e.g. Chrome 

##### Analyzing taxonomic composition of the datasets with phyloflash: mapping of reads to small subunit ribosomal RNA database 
## and targeted assembly SSU rRNA sequences #####

conda deactivate
conda activate communityComposition

phyloFlash.pl -man

## Defining the location of the SSU rRNA database 
##(SILVA non-redundnat reference sequences, edition 138.1, pre-processed and formatted for phyloFlash)
phyloDB=/vol/funmic/databases/phyloflashSilvaDB
READ_DIRECTORY=/vol/funmic/datasets/kelpBiofilm
OUTPUTDIRECTORY=/vol/funmic/phyloFlashOut/

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
			-CPUs 12 > /vol/funmicoflashlog.txt
		
		mv PhyloFlash_Kelp_${SAMPLE}* ${OUTPUTDIRECTORY}/${SAMPLE}
	done
	
##Time: 45 min ##
conda deactivate

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
	
## Time with 7 cores: 53 min for k21, 30 min for subsequent iterations, 9h 40 min per metagenome


## Simplify contig names down to simple IDs, removing everything that follows the space in the header
sed -i 's/\ .*//g' ${OUTPUTDIRECTORY}/*/*.contigs.fa
for ASSEMBLY in $(ls ${OUTPUTDIRECTORY} | xargs -n 1 basename);
	do
		sed -i 's/>/>'$ASSEMBLY'_/g' ${OUTPUTDIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa
	done

## okay, let's look at it with quast
conda deactivate

conda activate assemblyQC

mkdir /vol/funmic/Kelp/assemblyQC
ASSEMBLY_DIRECTORY=/vol/funmic/Kelp/assemblies/
OUTPUTDIRECTORY=/vol/funmic/Kelp/assemblyQC/

for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY}/ | xargs -n 1 basename);
	do 
		quast -t 12 \
			-o ${OUTPUTDIRECTORY}/${ASSEMBLY}/ \
			${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
			> /vol/funmic/Kelp/logs/${ASSEMBLY}_quastLog.txt
	done

conda deactivate	
## trivial amount of time and memory used.

####### Read mapping for differential coverage #########
conda activate binning

mkdir /vol/funmic/Kelp/mapping
ASSEMBLY_DIRECTORY=/vol/funmic/Kelp/assemblies/
READ_DIRECTORY=/vol/funmic/datasets/kelpBiofilm
OUTPUTDIRECTORY=/vol/funmic/Kelp/mapping

## each of the three read sets is mapped to each of the three assemblies: 9 mapping operations in total
## minid specifies the approximate nucleotide sequence identity to aim for when mapping reads. 
## Real identity of the alignment may be slightly higher or lower.
## idfilter parameter sets the real cut-off for reported alignments/mappings. Reads that map with an nucleotide identity less than 95% won't be reported.
## 95% is the average nucleotide identity between microbial genomes of the same species (e.g. different strains of the same species).

for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY} | xargs -n 1 basename);
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
				> /vol/funmic/Kelp/logs/${SAMPLE}_to_${ASSEMBLY}_bbmaplog.txt
				
				bash ${OUTPUTDIRECTORY}/${SAMPLE}_to_${ASSEMBLY}.sh
				rm ${OUTPUTDIRECTORY}/*.sam
			done
	done

##Time: 3h 20 min

####### BINNING ######



mkdir /vol/funmic/Kelp/binning

#defining some variables. DATA_FOLDER should be an absolute path to where the data for the project is stored, with fastq reads being in a subfolder called "reads".
#PROJECTNAME is a prefix which was used for the assembly etc.
#METAWATTPATH should be the absolute path to the folder with the Metawatt jar file. METAWATT_DATABASES is the absolute path to the Metawatt databases
#CONCOCTPATH should be the absolute path to the concoct folder containing subfolders "bin" and "scripts"

MAPPING_DIRECTORY=/vol/funmic/Kelp/mapping
ASSEMBLY_DIRECTORY=/vol/funmic/Kelp/assemblies
BINNING_DIRECTORY=/vol/funmic/Kelp/binning


#Binning of the metagenomes co-assembly 

#This script is part of the MetaBAT module and reads out the bam files to generate a coverage table (incl. coverage variation) for the contigs

for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY}/ | xargs -n 1 basename);
	do  
		mkdir ${BINNING_DIRECTORY}/${ASSEMBLY}
		mkdir ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT
		jgi_summarize_bam_contig_depths \
			--outputDepth ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_coverage_Depths.txt \
			--referenceFasta  ${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
			${MAPPING_DIRECTORY}/*_to_${ASSEMBLY}_sorted.bam
    done
    
#run MetaBAT
for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY}/ | xargs -n 1 basename);
	do     
	metabat2 \
        -i ${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
        -a ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_coverage_Depths.txt \
        -o ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_MetaBAT \
        -t 12
	done

## Time neglidgible (few minutes)

## perform binning with MaxBin
for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY}/ | xargs -n 1 basename);
  do 
    mkdir ${BINNING_DIRECTORY}/${ASSEMBLY}/MaxBin
# create MaxBin input tables by extracting the coverage columns from the metabat input table and saving as separate files. Header line is removed in the process
    for i in $(seq 4 +2 8); 
        do
            cut -f 1,"$i"  ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_coverage_Depths.txt | perl -0pe 's/contigName.*\n//g' > ${BINNING_DIRECTORY}/${ASSEMBLY}/MaxBin/${ASSEMBLY}_abundance_"$i".txt
        done
  
    ls ${BINNING_DIRECTORY}/${ASSEMBLY}/MaxBin/*.txt > ${BINNING_DIRECTORY}/${ASSEMBLY}/MaxBin/${ASSEMBLY}_abundances_list.txt
done

# run MaxBin
for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY}/ | xargs -n 1 basename);
  do    
    run_MaxBin.pl \
        -contig ${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
        -abund_list ${BINNING_DIRECTORY}/${ASSEMBLY}/MaxBin/${ASSEMBLY}_abundances_list.txt \
        -out ${BINNING_DIRECTORY}/${ASSEMBLY}/MaxBin/${ASSEMBLY}_MaxBin \
        -thread 14
  done


## Time: 2 h 7 min.

## perform binning with CONCOCT

MAPPING_DIRECTORY=/vol/funmic/Kelp/mapping
ASSEMBLY_DIRECTORY=/vol/funmic/Kelp/assemblies/
BINNING_DIRECTORY=/vol/funmic/Kelp/binning

## create contig coverage tables for CONCOCT from the MetaBAT coverage tables.
## The first column with coverage is actually column 4.
## Start from column 4, pick out every second column (CONCOCT doesn't need the coverage variation columns, only the coverage itself) until you got all the samples. 
## Also remove the first line with the headers.

for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY} | xargs -n 1 basename);
	do 
		mkdir ${BINNING_DIRECTORY}/${ASSEMBLY}/concoct
		cut -f 1,4,6,8 ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_coverage_Depths.txt > ${BINNING_DIRECTORY}/${ASSEMBLY}/concoct/${ASSEMBLY}_coverage_Depths.txt
	done

for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY} | xargs -n 1 basename);
	do     
	concoct \
		-t 14 \
        --composition_file ${ASSEMBLY_DIRECTORY}/${ASSEMBLY}/${ASSEMBLY}.contigs.fa \
        --coverage_file ${BINNING_DIRECTORY}/${ASSEMBLY}/concoct/${ASSEMBLY}_coverage_Depths.txt \
		-r 250 \
        -b ${BINNING_DIRECTORY}/${ASSEMBLY}/concoct/${ASSEMBLY}_CONCOCT \
        --no_original_data
	done

## Time: 1:30 h

## Choosing best bins from the results of different binning tools for each assembly with DAStool. DAStool judges bin size, completeness and contamination values
## generate DAS_Tool input files: Tables with two columns assigning each contig to a putative genomic bin. For each assembly we need such tables from every binning tool.

## Attribution of contigs to bins according to MetaBAT
for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY}/ | xargs -n 1 basename);
  do 
	mkdir ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool
	  for BIN in $(ls ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/${ASSEMBLY}_MetaBAT.*.fa | xargs -n 1 basename | sed 's/\.fa//g'); 
		do 
			grep '^>' ${BINNING_DIRECTORY}/${ASSEMBLY}/MetaBAT/"$BIN".fa | sed "s/\t.*/\t$BIN/g" | sed 's/>//g' \
					>> ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool/${ASSEMBLY}_MetaBAT_binning.txt 
		done
  done

## Attribution of contigs to bins according to MaxBin
for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY} | xargs -n 1 basename);
  do 
	  for BIN in $(ls ${BINNING_DIRECTORY}/${ASSEMBLY}/MaxBin/*.fasta | xargs -n 1 basename | sed 's/\.fasta//g'); 
		do 
			grep '^>' ${BINNING_DIRECTORY}/${ASSEMBLY}/MaxBin/${BIN}.fasta | sed "s/$/\t$BIN/g" | sed 's/>//g' \
					>> ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool/${ASSEMBLY}_MaxBin_binning.txt 
		done
  done


## Attribution of contigs to bins according to CONCOCT
for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY} | xargs -n 1 basename);
  do 
	sed '1d' ${BINNING_DIRECTORY}/${ASSEMBLY}/concoct/${ASSEMBLY}_CONCOCT_clustering_gt1000.csv | \
	perl -0pe "s/\,/\t$ASSEMBLY\_CONCOCT/g" \
	> ${BINNING_DIRECTORY}/${ASSEMBLY}/DAStool/${ASSEMBLY}_CONCOCT_binning.txt
  done

conda deactivate

## Run DAS_Tool to generate consensus bins from the results of different binners
conda activate refine

ASSEMBLY_DIRECTORY=/vol/funmic/Kelp/assemblies
BINNING_DIRECTORY=/vol/funmic/Kelp/binning

for ASSEMBLY in $(ls ${ASSEMBLY_DIRECTORY} | xargs -n 1 basename);
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

conda activate dereplication

mkdir /vol/funmic/Kelp/binning/drep
mkdir /vol/funmic/Kelp/binning/drep/all_bins
mkdir /vol/funmic/Kelp/binning/drep/dereplicated_bins

BINNING_DIRECTORY=/vol/funmic/Kelp/binning
ALL_BINS_DIRECTORY=/vol/funmic/Kelp/binning/drep/all_bins
DEREP_BINS_DIRECTORY=/vol/funmic/Kelp/binning/drep/dereplicated_bins

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

mkdir /vol/funmic/Kelp/Final_bins

FINAL_BINS_DIRECTORY=/vol/funmic/Kelp/Final_bins

mkdir ${FINAL_BINS_DIRECTORY}/contigs
scp ${DEREP_BINS_DIRECTORY}/*.fa ${FINAL_BINS_DIRECTORY}/contigs/

## Time 1 h 20 min

conda deactivate

## Check MAG completeness and contamination with the new checkm2. Results from CheckM 1.2.2 should be available from dRep

conda activate completeness_estimate

BINNING_DIRECTORY=/vol/funmic/Kelp/binning
# ALL_BINS_DIRECTORY=/vol/funmic/Kelp/binning/drep/all_bins

DEREP_BINS_DIRECTORY=/vol/funmic/Kelp/binning/drep/dereplicated_bins/dereplicated_genomes

checkm2 predict \
	--allmodels \
	--threads 12 \
	--remove_intermediates \
	--input ${FINAL_BINS_DIRECTORY}/contigs \
	--extension .fa \
	--output-directory ${FINAL_BINS_DIRECTORY}/CheckM
	
## Time 6 min
conda deactivate

##  Classify the bins taxonomically with GTDB-Tk
conda activate assignTaxonomy

BINNING_DIRECTORY=/vol/funmic/Kelp/binning
DEREP_BINS_DIRECTORY=/vol/funmic/Kelp/binning/drep/dereplicated_bins/dereplicated_genomes


gtdbtk classify_wf \
    --mash_db /vol/funmic/databases/gtdb/release220/gtdb-tk-r220.msh \
	--cpus 12 \
	--pplacer_cpus 12 \
	--genome_dir {FINAL_BINS_DIRECTORY}/contigs/ \
	-x fa \
	--out_dir {FINAL_BINS_DIRECTORY}/GTDB
	
## Time 4h 45 min

##-------------------------------------------------------------------------------------------------------------------------
## Functional prediction with MicroTrait
##-------------------------------------------------------------------------------------------------------------------------

conda activate annotation

cd ${FINAL_BINS_DIRECTORY}

/usr/bin/R

library(microtrait)
require(stringr)

genomes_dir <- "contigs/"
out_dir <- "microtrait/"

genomes <- str_replace_all(list.files(genomes_dir, pattern=".fa"), ".fa$", "")

for(GENOME in genomes) {
  results_protein <- extract.traits(in_file=paste(genomes_dir, GENOME, ".fa", sep=""),
                                    out_dir = "microtrait/",
                                    type="genomic",
                                    mode="single",
                                    growthrate_predict = T,
                                    optimalT_predict = T)
  
    write.table(results_protein$microtrait_result$trait_counts_atgranularity3,
              file=paste(out_dir,GENOME,"trait_counts_atgranularity3.tsv", sep=""),
              sep="\t",
              quote=F,
              row.names=F)
}


rds_files <- list.files(out_dir, 
                        pattern=".rds", 
                        full.names=T)

genomeset_results = make.genomeset.results(rds_files,
                                           ids = genomes,
                                           growthrate = T,
                                           optimumT = T,
                                           ncores = 12)


write.table(genomeset_results$hmm_matrix, 
                file=paste(out_dir,"Genomes_hmm_matrix.tsv", sep=""), 
                sep="\t", 
                row.names=F, 
                quote=F, 
                na="")                                           

write.table(genomeset_results$trait_matrixatgranularity3, 
                file=paste(out_dir,"Genomes_trait_matrixatgranularity3.tsv", sep=""), 
                sep="\t", 
                row.names=F, 
                quote=F, 
                na="")

write.table(genomeset_results$rule_matrix, 
                file=paste(out_dir,"Genomes_rule_matrix.tsv", sep=""), 
                sep="\t", 
                row.names=F, 
                quote=F, 
                na="")


quit(save="no")


# Annotation of proteins with InterPro

conda activate annotation

FINAL_BINS_DIRECTORY=/vol/funmic/Kelp/Final_bins

cd ${FINAL_BINS_DIRECTORY}

mkdir proteins

for MAG in $(ls contigs/*.fa | xargs -n 1 basename | sed 's/.fa//g');

do
     prodigal \
        -i contigs/${MAG}.fa \
        -a proteins/${MAG}.faa \
        -q
done

# prodigal generates amino acid sequences which end with an aesterix, *, symbolizing the stop codon.
# InterProScan does not accept such seqeunces. The aesterix has to be removed with search & replace

sed -i 's/\*//g' proteins/*.faa 

mkdir Annotation

INTERPRO_DATA_DIR=/vol/funmic/databases/interproscan-5.73-104.0

for MAG in $(ls proteins/*.faa | xargs -n 1 basename | sed 's/.faa//g');

do

    ${INTERPRO_DATA_DIR}/interproscan.sh \
		-cpu 12 \
        -appl Pfam,TIGRFAM,SUPERFAMILY \
		-i proteins/${MAG}.faa \
		-f tsv \
		-o Annotation/${MAG}_InterProScan.tsv

done
       
# Annotation of proteins with EggNOG-Mapper (VERY SLOW)


for MAG in $(ls proteins/*.faa | xargs -n 1 basename | sed 's/.faa//g');
do

    emapper.py \
        --cpu 12 \
        -m diamond \
        --qtype seq \
        --seed_ortholog_evalue 0.00001 \
        -i proteins/${MAG}.faa \
        --output_dir Annotation/ \
        -o ${MAG}_EggNOG
    
done

  sed -i '/^##/d' .*_EggNOG.emapper.annotations







####-----------------------------------------------------------------------------------------------------------------------
#### Generating an Anvi'O database (right now Anvi'O interactive mode is not working)
####-----------------------------------------------------------------------------------------------------------------------

mkdir /vol/funmic/Anvio_kelp
mkdir /vol/funmic/Anvio_kelp/input

ANVIO_DIRECTORY=/vol/funmic/Anvio_kelp
DEREP_BINS_DIRECTORY=/vol/funmic/binning/drep/dereplicated_bins/dereplicated_genomes
READ_DIRECTORY=/vol/funmic/datasets/kelpBiofilm

## Combine all bins into one FASTA FILE

cat ${DEREP_BINS_DIRECTORY}/*.fa > ${ANVIO_DIRECTORY}/input/All_bins.fasta

## Map the reads from all samples to the combined fasta file

conda activate RawReadProcessing

ANVIO_DIRECTORY=/vol/funmic/Anvio_kelp
DEREP_BINS_DIRECTORY=/vol/funmic/binning/drep/dereplicated_bins/dereplicated_genomes
READ_DIRECTORY=/vol/funmic/datasets/kelpBiofilm

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

ANVIO_DIRECTORY=/vol/funmic/Anvio_kelp

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

DEREP_BINS_DIRECTORY=/vol/funmic/binning/drep/dereplicated_bins/dereplicated_genomes

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

Kaiju_DB=/vol/funmic/anaconda3/envs/anvio-8/data

  kaiju \
    -z 24 \
    -E 0.00001 \
    -e 5 \
    -t ${Kaiju_DB}/nodes.dmp \
    -f ${Kaiju_DB}/kaiju_db_refseq_nr.fmi \
    -i ${ANVIO_DIRECTORY}/Kelp.faa \ 
	-p \
    -o ${ANVIO_DIRECTORY}/Input/"$SAMPLE".kaiju


		
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

