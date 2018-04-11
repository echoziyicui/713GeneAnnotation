#Author: Saideep Gona

import sys
import os
import subprocess
import multiprocessing

# Inputs
    # transcriptome_path (not an argument, assume the file is in the step1.1 directory from previous step)

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

run_single_augustus(transcriptome_path, 'human', 'augustus_gff.gff')




