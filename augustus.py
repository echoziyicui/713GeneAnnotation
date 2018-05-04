#Author: Saideep Gona

import sys
import os
import subprocess
import multiprocessing

# Dependencies
# blat
# augustus

# Inputs
    # transcriptome_path (not an argument, assume the file is in the step1.1 directory from previous step)

# Outputs
    # augustus_output.gff

cwd = os.getcwd()
transcriptome_path = cwd + "/step1.1/intermediate_fasta.fasta"
reference_directory = cwd + "/references/"
intermediate_dir = cwd + "/step1.2/"

try:
    os.mkdir(intermediate_dir)
except:
    pass

to_annotate = sys.argv[1]
species = sys.argv[2]
# ref_path = reference_directory + reference

# blat -minIdentity=92 genome.fa cdna.fa cdna.psl
# blat2hints.pl --in=/step1.2/cdna.psl --out=step1.2/hints.E.gff
# augustus --species=human --hintsfile=hints.E.gff --extrinsicCfgFile=extrinsic.ME.cfg genome.fa

blat_1 = "blat -minIdentity=92 '"+to_annotate+"' '"+transcriptome_path+"' '"+intermediate_dir+"transcriptome.psl'"
blat_2 = "blat2hints.pl --in='"+intermediate_dir+"transcriptome.psl'"+" "+"--out='"+intermediate_dir+"hints.gff'"

def run_single_augustus(transcriptome_path, species, output_file):

    augustus_args = [
        "augustus",
        "--species=" + species,
        "--hintsfile=" +intermediate_dir+"hints.gff",
        "--extrinsicCfgFile=extrinsic.ME.cfg",
        to_annotate
    ]

    out = open(intermediate_dir + output_file,"w")

    augustus = subprocess.Popen(augustus_args, stdout=out)

    return
print("Starting Step 1.2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print("Running Blat Transcriptome Aligner, Outputting .psl")
os.system(blat_1)
print("Generating Rnaseq Hints with blat2hints.pl, Outputting hints.gff")
os.system(blat_2)
print("Annotating Genome with Augustus, Outputting .gff")
run_single_augustus(transcriptome_path, species, 'augustus_gff.gff')
print("Step 1.2 Complete >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")




