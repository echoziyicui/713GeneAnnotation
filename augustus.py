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

cwd = os.getcwd()
transcriptome_path = cwd + "/step1.1/intermediate_fasta.fasta"

intermediate_dir = cwd + "/step1.2/"

try:
    os.mkdir(intermediate_dir)
except:
    pass

def run_single_augustus(transcriptome_path, species, output_file):

    augustus_args = [
        "augustus",
        "--species=" + species,
        transcriptome_path
    ]

    out = open(intermediate_dir + output_file,"w")

    augustus = subprocess.Popen(augustus_args, stdout=out)

    return 

run_single_augustus(transcriptome_path, 'human', 'intermediate_gff.gff')




