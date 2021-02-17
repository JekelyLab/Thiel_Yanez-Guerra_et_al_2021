This bash script was created on a Linux Ubuntu 18.04 operating system.

The script performs a tBLASTn search with any number of genes on any number of transcriptomes.
The candidate nucleotide sequences get translated into amino acid sequences and those that can't get translated (e.g. due to an indel) will be stored separately.

Add this script in an extra folder and create the following folders inside that folder: "Apps", "Transcriptomes", "Reference_genes".

Collect transcriptomes, change name of transcriptomes to Speciesname.fasta. It is important that the files end in ".fasta". Move transcriptome files into folder "Transcriptomes". The file name will become part of the sequence headers. (Keep a copy of the transcriptomes in an extra folder; this script changes the original headers and adds the file name everytime it is run on a transcriptome.) 
Create files with reference precursors (amino acid sequences); each file can have multiple query sequences of the same neuropeptide; name files according to the neuropeptide and ending with ".fasta"; move query files into "Reference_genes". (Alternatively you can also use a single file with different peptides, but then they have to be sorted manually.)


Dependencies:

- Perl
- BLAST commandline tool from ncbi (in PATH)
- Transdecoder from https://github.com/TransDecoder (add into "Apps/TransDecoder/")
- fastagrep.pl script from https://github.com/rec3141/rec-genome-tools/blob/master/bin/fastagrep.pl (add into "Apps/")


Make fastagrep.pl executable
Make this script executable
Execute this script from the command line once everything is set up as described


Candidate sequences are in the folder "Output/".
Sequences that were detected by tBLASTn but could not get translated are in the corresponding gene folder in "Output/protein/" named 003_[gene_name]_all_missing_ORFs_nuc_seqs.txt'"



