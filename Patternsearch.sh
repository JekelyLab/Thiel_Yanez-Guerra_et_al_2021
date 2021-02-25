#!/bin/bash

# This bash script was created on a Linux Ubuntu 18.04 operating system.
# This script uses bash regular expressions to search for repetitive cleavage and amidation sites in amino acid sequences.

# Create a folder called "Transcriptomes" in the same folder as this script.
# Add the translated Transcriptomes (amino acids sequences) for scanning into the folder "Transcriptomes".
# It is important that the transcriptome files end in ".fasta"

# Make this script executable

# NOTE! There will be a lot of false positives if this script is run on a whole transcriptome. 
	# It is therefore resommended to first sort for only those sequences that have a signal peptide with another program such as SignalP (https://services.healthtech.dtu.dk/service.php?SignalP-5.0).


mkdir "./Temp"

# the first loop creates single line fasta files
for sequencecollection in Transcriptomes/*
do
	Filename="${sequencecollection#Transcriptomes/}"

	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < "${sequencecollection}" > "Temp/${Filename}"
done

# the second loop uses the new single line fasta files to scan sequences for the specified repetitions
for transcriptome in Temp/*
do
	transcriptome_file="${transcriptome#Temp/}"
	Species="${transcriptome_file%.fasta}"
	MCP_candidates="${Species}_MCP_candidates.fasta"

	grep -E -B 1 -e "((.{2,25}G[KR][KR])|(.{2,25}[KR](.{1}|.{3}|.{5})G[KR*])){3,}" -e "(.{5,25}[KR][KR]){5,}" "${transcriptome}" | grep -E -v "[-][-]" > "${MCP_candidates}"

done

rm -R "./Temp"

# Explanation of regular expressions

	# -E option = interpret patterns as extended regular experssions
	# -B 1 option = print additionally 1 line from above the context (in this case the sequence header)
	# -e option = when searching for multiple Patterns
	# [] square brackets contain the list of characters to scan for: [KR] means either K or R
	# . dot = any character
	# .* dot star = any number (including 0) of any character
	# () brackets indicate a group that belongs together
	# number in {} curly brackets indicate the range of how many of the previous character
	# | pipe redirects the output and uses the following command on it
	# -v option = check for lines that do NOT have the following expression (in this case it skips the -- that is introduced in the output for some reason
	# if special characters like { or * should be included in the pattern, then they should be put into square brackets []

# Further explanations at:
	# https://www.cyberciti.biz/faq/grep-regular-expressions/
	# https://linuxtechlab.com/bash-scripting-learn-use-regex-basics/
	# https://www.gnu.org/savannah-checkouts/gnu/grep/manual/grep.html#Command_002dline-Options


# alternative sequence to scan for:
	# -e "(([KR].{2,25}G[KR][KR])|([KR].{2,25}[KR](.{1}|.{3}|.{5})G[KR*])){3,}"
	# or
	# "G[KR][KR].{2,25}G[KR][KR].{2,25}G[KR][KR]"
	
# or simpler
	# grep -E -B 1 -e "G[KR][KR].{2,50}G[KR][KR].{2,50}G[KR][KR]" -e "[KR](.{2}|.{4}|.{6})R.{2,25}G[KR][KR]" -e "(([KR].{2,25}G[KR][KR])|([KR].{2,25}[KR](.{2}|.{4}|.{6})G[KR])|.{3,25}[KR][KR]){2,}" "${transcriptome}" | grep -E -v "[-][-]" > output.fasta
