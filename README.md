# 713GeneAnnotation
Group project for 03-713 <br>
This project is an RNA-seq-based gene annotation, verified by ribosome-profiling-based translatome identification.
Instead of using the traditional methods for gene annotation, for example, RNA-seq profiling or ribosome profilling
This pipeline is a combination of using a gene annotation software called Augustus, which states as a pre-trained hidden markov modle to increase the produced gene sequence annotation's accuracy. 
By transforming the whole genome sequence into the Augustus, this pipeline, different from others, will be able to produce it's own GTF file, which stands for a file has the locations corresponding to the specific gene transcription features in the whole DNA sequence. 
The main purpose of incorporation Augustus into this pipeline is to increase our accuracy of gene annotation predictions with a given RNA-seq data as a useful "hints", since the RNA-seq is able to provide clues about underlying gene structure, which could contributes to a higher confidence in our prediction. 
By reprocessed our gene annotation sequence through Augustus, we finally will be able to have the GTF file that is the input for our ribo-code software. And by using the ribosome profiling data, we are more likely to explore the novel gene annotation location.
