#Author: Saideep Gona

import sys
import os
import subprocess
import multiprocessing

# Inputs
    # reference_genome: filename of reference genome
    # input_fastqs: filenames of all input fastqs

# Outputs
    # augustus_output.gff 


reference_genome = sys.argv[1]                 # Name of reference genome - eg. hg38.fa - 
input_fastqs = sys.argv[2:-1]               # Names of input fastqs
cores = multiprocessing.cpu_count()

print("Number of cpus: ", cores)

cwd = os.getcwd()
reference_directory = cwd + "/references/" 
fastq_directory = cwd + "/fastq/"
ref_genome = reference_directory + reference_genome

intermediate_dir = cwd + "/step1/"
try:
    os.mkdir(intermediate_dir)
except:
    pass


# Arguments for bwa alignment call
bwa_args = [                                
    "bwa",
    "mem",
    "-t",
    cores,
    ref_genome
] + input_fastqs

# Arguments for samtools sort call (piped from bwa mem)
sam_sort_args = [
    "samtools",
    "sort",
    "-O",
    "BAM",
    "-o",
    intermediate_dir + "intermediate_bam.bam"
]

# Arguments for scallop assembly call
scallop_args = [
    "scallop",
    "-i",
    "intermediate_bam.bam",
    "-o",
    intermediate_dir + "intermediate_gtf.gtf"
] 

# Arguments for gffread gtf conversion to fasta call
gffread_args = [
    "gffread",
    "-g",
    ref_genome,
    cwd + "intermediate_gtf.gtf",
    "-w",
    intermediate_dir + "intermediate_fasta.fasta"
]




print("Starting Pipeline")
print("Mapping Reads to Reference")

bwa_error = open("bwa_error.txt", 'w')
bwa = subprocess.Popen(bwa_args, stdout=subprocess.PIPE, stderr=bwa_error)

sam_sort_error = open("sam_sort_error.txt", 'w')
sam_sort = subprocess.Popen(sam_sort_args, stdin=bwa.stdout, stderr=bwa_error)
bwa.stdout.close()
produce_bam = sam_sort.communicate()
bwa.wait()

scallop = subprocess.Popen(scallop_args)

gffread_args = subprocess.Popen(gffread_args)