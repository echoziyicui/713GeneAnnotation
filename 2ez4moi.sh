genomefile=$1
inputfastq=$2
python3 raw_rna_seq_to_augustus.py $genomefile $inputfastq
python3 augustus.py
more step1.2/augustus_gff.gff |grep chr > step1.2/augustus_final_output.gtf