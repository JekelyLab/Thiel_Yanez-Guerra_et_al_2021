This bash script was created on a Linux Ubuntu 18.04 operating system.

This script uses bash regular expressions to search for repetitive cleavage and amidation sites in amino acid sequences.


Create a folder called "Transcriptomes" in the same folder as this script.
Add the translated Transcriptomes (amino acids sequences) for scanning into the folder "Transcriptomes".
It is important that the transcriptome files end in ".fasta"

Make this script executable

Note:
There will be a lot of false positives if this script is run on a whole transcriptome. It is therefore worth it to first sort for only those sequences that have a signal peptide with another program such as SignalP (https://services.healthtech.dtu.dk/service.php?SignalP-5.0)