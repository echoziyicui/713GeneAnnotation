#!/bin/bash
################################################################################
# A bash script for automating the process of RiboCode to annotate genome.
#
# Author:    Group D
# Version:   2.0
# Date:      04/26/2018
#
# Function:
# ---------
# 1. Remove ribosomal RNA derived read use Bowtie and STAR
# - Bowtie finds the non-rRNA reads in ribo-seq
# 2. Align the non-rRNA reads using STAR
# - STAR builds index for alignment
# - STAR aligns the non-rRNA ribo reads back to the reference genome
# 3. Identify translated ORFS in genome based on ribo-seq
# ---------
# Usage:
# ---------
# $1:    Genome fasta
# $2:    Genome GTF
# $3:    rRNA fasta
# $4:    ribo-seq fastq (ribosome profiling data)
# $5:    folder name to contain the star index
# ---------
# Test example used: mouse liver m16(GRCm38.p5)
# Command line:
# bash ribocode.sh GRCm38.p5.genome.fa gencode.vM16.annotation.gtf \
# mouse_pre_rRNA.fa mouse_rRNA.fastq folder
# ---------
# $1: <GRCm38.p5.genome.fa>
#      https://www.gencodegenes.org/mouse_releases/16.html
# $2: <gencode.vM16.annotation.gtf>
#     (GRCm38.p5, comprehensive annotation only for chromosome)
#      https://www.gencodegenes.org/mouse_releases/16.html
# $3: <mouse_pre_rRNA.fa> (45s rRNA)
#      http://rnacentral.org/rna/URS000075CEB1/10090
# $4: <mouse_rRNA.fastq> (mouse liver)
#      pre-downloaded using fastq-dump command of SRA_toolkit
#      https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1630812
# $5: <folder>
#      required by STAR for storing the index built
################################################################################
### File Preparation
# Re-name the command line arguments
genome_fasta=$1
genome_GTF=$2
rRNA_fasta=$3
ribo_fastq=$4
folder=$5
# Create an empty folder to hold the STAR index.
mkdir $folder
# Rename the ribo-seq to match the bowtie-index prefix.
mv $ribo_fastq mouse_rRNA.fastq
################################################################################
### Preprocessing
##Step1. Use bowtie to remove rRNA in the ribo-seq.
#export PATH=${PATH}:/usr1/home/sgona/bowtie-0.12.7

$bowtie_path/bowtie-build $rRNA_fasta mouse_rRNA
$bowtie_path/bowtie -p 8 -norc --un un_aligned_ribo_mouse.fastq -q mouse_rRNA \
mouse_rRNA.fastq mouse_rRNA.align

##Step2: Use star to align the clean reads to reference genome fasta file to
##       get a ribo-seq bam file.
# a. Build index for star.
#export PATH=${PATH}:/usr1/home/sgona/STAR/bin/Linux_x86_64_static
$star_path/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $folder \
--genomeFastaFiles $genome_fasta --sjdbGTFfile $genome_GTF
# b. Align the clean read back to fasta genome.
$star_path/STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 \
--genomeDir $folder --readFilesIn un_aligned_ribo_mouse.fastq  \
--outFileNamePrefix mouseliver --outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 \
--outFilterMatchNmin 16 --alignEndsType EndToEnd
################################################################################
### ORF Identification
#Run ribocode to find ORFs in the given genome fasta based on rRNA-removed
#ribo-seq bam file.
#export PATH=${PATH}:/usr1/home/sgona/RiboCode/RiboCode
$ribocode_path/RiboCode_onestep -g $genome_GTF -f $genome_fasta \
-r mouseliverAligned.toTranscriptome.out.bam -l no -o RiboCode_ORFs_result
################################################################################
### Postprosessing
# Re-format the Ribocode results (txt) to GTF form
awk 'BEGIN {OFS="\t"};{print $8,$2,$7,$13,$14,$28,$9,".","gene_id \""$5"\"; gene_name \""$6"\"; "}'  RiboCode_ORFs_result_collapsed.txt >ribocode_gtf.gtf


### ~/.profile ### in manual makesure they have the environment variable
