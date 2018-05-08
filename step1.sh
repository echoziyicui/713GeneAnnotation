genomename=$1
gtfname=$2
inputfastq=$3
to_annotate=$4
species=$5

python3 raw_rna_seq_to_augustus.py $genomename $gtfname $inputfastq
python3 augustus.py $to_annotate $species
more step1.2/augustus_gff.gff |grep chr > step1.2/augustus_final_output.gtf

awk 'BEGIN {OFS="\t"};{print $8,$2,$7,$13,$14,$28,$9,".","gene_id \""$5"\"; gene_name \""$6"\"; "}'  RiBoCode_ORFs_result.txt >ribocode_gtf.gtf
perl ribocode_gtf2gtf.pl
