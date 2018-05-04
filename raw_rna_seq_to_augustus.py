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
reference_gtf = sys.argv[2]                    # Name of reference gtf file
input_fastqs = sys.argv[3:]                    # Names of input fastqs
cores = str(multiprocessing.cpu_count()-1)


print("Number of cpus: ", cores)

cwd = os.getcwd()
reference_directory = cwd + "/references/"
fastq_directory = cwd + "/fastq/"
# if not os.path.isdir(cwd + fastq_directory):
#     print("no fastq directory")
#     sys.exit()

ref_genome = reference_directory + reference_genome
ref_gtf = reference_directory + reference_gtf

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

cat_str = "cat " + " ".join(input_fastqs) + " > " + fastq_directory+"big.fastq"

#star_str = '''STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir mouse_liver_star_index
# --readFilesIn '''+fastq_directory+'''big.fastq --outFileNamePrefix'''+ intermediate_dir+'''temp --outSAMtype BAM SortedByCoordinate --
# quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd'''
print(ref_genome)
star_ind = "STAR --runThreadN=8 --runMode=genomeGenerate --genomeDir="+reference_directory+ " --genomeFastaFiles="+ref_genome+ " --sjdbGTFfile="+ref_gtf

star_str1 = '''STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir '''
star_str2 = reference_directory+''' --readFilesIn '''+fastq_directory+'''big.fastq --outFileNamePrefix '''
star_str3 = intermediate_dir+'''temp --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM '''
star_str4 = '''GeneCounts --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd'''
star_str = star_str1+star_str2+star_str3+star_str4

rename_str = "mv " + intermediate_dir + "tempAligned.sortedByCoord.out.bam " + intermediate_dir + "intermediate_bam.bam"

#alignment_full = "bwa mem -t "+cores+ref_genome+" ".join(input_fastqs)+"|samtools sort -O BAM -o "+intermediate_dir+"intermediate_bam.bam"
#print(alignment_full)
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
print("Indexing Reference")
if not os.path.exists(ref_genome+".fai"):
    os.system(star_ind)

print("Mapping Reads to Reference, Outputting .bam")

os.system(cat_str)

os.system(star_str)

os.system(rename_str)

print('Constructing Transcriptome From Mapped Reads, Outputting .gff')

scallop = subprocess.run(scallop_args)

print('Constructing .fasta From Transcriptome')

gffread_args = subprocess.run(gffread_args)

print('Step 1.1 Complete >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

 #Modify intermediate_fasta.fasta to have chromosome names

