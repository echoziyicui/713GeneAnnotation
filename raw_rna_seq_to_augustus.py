#Author: Saideep Gona

import sys
import os
import subprocess
import multiprocessing
import glob

# Inputs
    # reference_genome: filename of reference genome
    # input_fastqs: filenames of all input fastqs

# Outputs
    # augustus_output.gff

# Dependencies
#  - BWA
#  - Samtools
#  - Scallop
#  - Gffread
#  - blat
#  - Augustus


reference_genome = sys.argv[1]                 # Name of reference genome - eg. hg38.fa -
input_fastqs = sys.argv[2:]               # Names of input fastqs
cores = str(multiprocessing.cpu_count()-1)


print("Number of cpus: ", cores)

cwd = os.getcwd()
reference_directory = cwd + "/references/"
fastq_directory = cwd + "/fastq/"
# if not os.path.isdir(cwd + fastq_directory):
#     print("no fastq directory")
#     sys.exit()

ref_genome = reference_directory + reference_genome

intermediate_dir = cwd + "/step1.1/"
try:
    os.mkdir(intermediate_dir)
except:
    pass

input_fastqs = [fastq_directory + fastq  for fastq in input_fastqs]
print(input_fastqs)
# input_fastqs = glob.glob(fastq_directory + "*.fastq")
# if len(input_fastqs) == 0:
#     print("no input fastqs")
#     sys.exit()

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

alignment_full = "bwa mem -t "+cores+ref_genome+" ".join(input_fastqs)+"|samtools sort -O BAM -o "+intermediate_dir+"intermediate_bam.bam"
print(alignment_full)
# Arguments for scallop assembly call
scallop_args = [
    "scallop",
    "-i",
    intermediate_dir + "intermediate_bam.bam",
    "-o",
    intermediate_dir + "intermediate_gtf.gtf"
]

# Arguments for gffread gtf conversion to fasta call
gffread_args = [
    "gffread",
    "-g",
    ref_genome,
    intermediate_dir + "intermediate_gtf.gtf",
    "-w",
    intermediate_dir + "intermediate_fasta.fasta"
]


print("Starting Pipeline ***************************************************************************")
print("Starting Step 1.1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print("Mapping Reads to Reference, Outputting .bam")

bwa_error = open(intermediate_dir + "bwa_error.txt", 'w')
print(bwa_args)
print(type(bwa_error))
bwa = subprocess.Popen(bwa_args, stdout=subprocess.PIPE, stderr=bwa_error)

sam_sort_error = open(intermediate_dir + "sam_sort_error.txt", 'w')
sam_sort = subprocess.Popen(sam_sort_args, stdin=bwa.stdout, stderr=bwa_error)
bwa.stdout.close()
produce_bam = sam_sort.communicate()
bwa.wait()

print('Constructing Transcriptome From Mapped Reads, Outputting .gff')

scallop = subprocess.run(scallop_args)

print('Constructing .fasta From Transcriptome')

gffread_args = subprocess.run(gffread_args)

print('Step 1.1 Complete >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

 #Modify intermediate_fasta.fasta to have chromosome names




