#!/bin/bash

# General: 
	# This bash script was created on a Linux Ubuntu 18.04 operating system.
	# The script performs a tBLASTn search with any number of genes on any number of transcriptomes.
	# Depending on the genes, the minimum e-values for positive hits will have to be adjusted at two places: in line 106 and line 157.
	# The candidate nucleotide sequences that were detected by tBLASTn get automatically translated into amino acid sequences and those that can't get translated (e.g. due to an indel) will be stored separately in "Output/protein/" named 003_[gene_name]_all_missing_ORFs_nuc_seqs.txt'"

# Dependencies:
	# - Perl
	# - BLAST commandline tool from ncbi (in PATH)
	# - Transdecoder from https://github.com/TransDecoder (add into "Apps/TransDecoder/")
	# - fastagrep.pl script from https://github.com/rec3141/rec-genome-tools/blob/master/bin/fastagrep.pl (add into "Apps/")

# To execute this script:
	# Add this script into an extra folder and create the following folders inside this folder: "Apps", "Transcriptomes", "Reference_genes".
	# Collect transcriptomes, change name of transcriptomes to Speciesname.fasta. It is important that the files end in ".fasta"
			# The file name will become part of the sequence headers. Keep a copy of the transcriptomes in an extra folder; this script changes the original headers and adds the file name everytime it is run on a transcriptome.
	# Move transcriptome files into folder "Transcriptomes".
	# Create separate files with reference genes (amino acid sequences); name files according to the gene with file extension ".fasta"
		# each file can have multiple query sequences of the same neuropeptide; alternatively you can also use a single file with different peptides, but then they have to be sorted manually.  
	# Move query files into "Reference_genes".

	# Make fastagrep.pl executable (e.g. chmod 777 fastagrep.pl)
	# Make this script executable (e.g. chmod 777 tBLASTn_script.sh)

	# Execute this script from the command line after everything is set up as described





# ---------- PART 1 ----------

	# the first part changes the sequence headers: it deletes all special characters, keeps the identifier and adds the file name. This is neccessary for fastagrep.pl script.
	# the first part also creates the BLAST database.

echo " "

for Transcriptome in Transcriptomes/*
do
	Transcriptome_file="${Transcriptome#Transcriptomes/}"
	Species="${Transcriptome_file%.fasta}"
	
	echo "Changing headers of ${Transcriptome_file}"
	
	sed -i "s, ,\t,g" "${Transcriptome}"
	awk '{FS = "\t"} {print $1}' ${Transcriptome} > tmp && mv tmp ${Transcriptome}
	sed -i "s,\t,__,g" "${Transcriptome}"
	sed -i -e 's.\+._.g' -e 's.\:._.g' -e 's.\;._.g' -e 's.\/._.g' -e 's.|._.g' -e 's.\\._.g' -e 's.{._.g' -e 's.\[._.g' -e 's.(._.g' -e 's.}._.g' -e 's.]._.g' -e 's.)._.g' -e 's.\^._.g' -e 's.\$._.g' -e 's.,._.g' -e 's.?._.g' -e 's.*..g' "${Transcriptome}"
	sed -i "s,>,>"${Species}"_,g" "${Transcriptome}"
	
	echo "Creating BLAST database for ${Species}"
	makeblastdb -dbtype nucl -in "${Transcriptome}" -out "./Transcriptomes/${Species}"
	echo " "
done

	mkdir "./Output"
	mkdir "./Output/protein"
	mkdir "./Output/nucleotide"
	mkdir "./Output/protein/tmp"
	mkdir "./Output/tmp"



# ---------- PART 2 ----------

	# the second part BLASTs each reference gene file against each transcriptome.
	# the outer loop goes through the reference gene files
	# the inner loop BLASTs the current reference gene file against all transcriptome files; 
		# the inner loop also translates the positive hits into amino acid sequences and stores those that did not get translated in a separate file (when a gene was detected by tBLASTn, but could not be translated into amino acids due to an indel)
	# the result is a separate file for each gene and each species, both as nucleotide as well as translated amino acid sequences, with all files from one gene saved in a common output folder.
		# there will also be a single file containing all genes of all transcriptomes together called "All_neuropeptide_prot_candidates.fasta"
	
for gene in Reference_genes/*
do
	Reference_file="${gene#Reference_genes/}"
	Reference_gene="${Reference_file%.fasta}"

	# create a separate folder for each gene in the output folder
	mkdir "./Output/nucleotide/${Reference_gene}"
	mkdir "./Output/nucleotide/${Reference_gene}/tBLASTn_output"
	mkdir "./Output/nucleotide/${Reference_gene}/tBLASTn_output/IDs"
	mkdir "./Output/protein/${Reference_gene}"
	mkdir "./Output/protein/${Reference_gene}/BLASTp_output"
	mkdir "./Output/protein/${Reference_gene}/BLASTp_output/IDs"
	mkdir "./Output/protein/${Reference_gene}/"
	mkdir "./Output/protein/${Reference_gene}/unclear_ORFs"
	mkdir "./Output/protein/${Reference_gene}/missing_ORFs"
	mkdir "./Output/protein/${Reference_gene}/missing_ORFs/tmp"

	
	for Transcriptome in Transcriptomes/*.fasta
	do
		Transcriptome_file="${Transcriptome#Transcriptomes/}"
		Species="${Transcriptome_file%.fasta}"
		tBLASTn_outputfile="${Species}_${Reference_gene}_tBLASTn_output.fasta"
		tBLASTn_output_IDs="${Species}_${Reference_gene}_tBLASTn_IDs.fasta"
		tBLASTn_candidates="${Species}_${Reference_gene}_nuc.fasta"
		
		echo "tBLASTn: ${Reference_gene} sequences against ${Species} transcriptome" 
		
		
# ---------- e-value adjustment 1 ----------		
	# this is the first place where the minimum e-value for the tBLASTn search has to be adjusted:
		tblastn -query "${gene}" -db "./Transcriptomes/${Species}" -task tblastn -evalue 1e-5 -out "./Output/nucleotide/${Reference_gene}/tBLASTn_output/${tBLASTn_outputfile}" -outfmt 6
		
		mkdir "./Output/nucleotide/${Reference_gene}/tmp1"
		cut -f 2 "./Output/nucleotide/${Reference_gene}/tBLASTn_output/${tBLASTn_outputfile}" > "./Output/nucleotide/${Reference_gene}/tmp1/${tBLASTn_outputfile}"

		mkdir "./Output/nucleotide/${Reference_gene}/tmp2"
		sort "./Output/nucleotide/${Reference_gene}/tmp1/${tBLASTn_outputfile}" > "./Output/nucleotide/${Reference_gene}/tmp2/${tBLASTn_outputfile}"
		rm -R "./Output/nucleotide/${Reference_gene}/tmp1"

		uniq "./Output/nucleotide/${Reference_gene}/tmp2/${tBLASTn_outputfile}" > "./Output/nucleotide/${Reference_gene}/tBLASTn_output/IDs/${tBLASTn_output_IDs}"
		rm -R "./Output/nucleotide/${Reference_gene}/tmp2"
		
		./Apps/fastagrep.pl -f "./Output/nucleotide/${Reference_gene}/tBLASTn_output/IDs/${tBLASTn_output_IDs}" -X "${Transcriptome}" > "./Output/nucleotide/${Reference_gene}/${tBLASTn_candidates}"
		sed -i -e "s,\.p,_y_y_y_y_p,g" "./Output/nucleotide/${Reference_gene}/${tBLASTn_candidates}"
		


		# the detected nucleotide sequences will be translated into protein sequences:
			# this is done by first translating the sequences into all ORFs and then performing a second BLAST (BLASTp) to detect the correct ORF.
			
		echo "Translating ${Species} ${Reference_gene} candidate nucleotide sequences into amino acid sequences"
		./Apps/TransDecoder/TransDecoder.LongOrfs -m 50 -t ./Output/nucleotide/${Reference_gene}/${tBLASTn_candidates}

		Candidate_ORFs="${tBLASTn_candidates%.fasta}_ORFs.fasta"
		Candidate_ORFs_db="${Candidate_ORFs%.fasta}"

		rm -R "${tBLASTn_candidates}.transdecoder_dir.__checkpoints_longorfs"
		rm *.cmds
		cp "${tBLASTn_candidates}.transdecoder_dir/longest_orfs.pep" "./Output/protein/tmp/${Candidate_ORFs}"
		sed -i -e "s, ,\t,g" -e "s,\.p,_x_x_x_x_p,g" "./Output/protein/tmp/${Candidate_ORFs}"
		awk '{FS = "\t"} {print $1}' ./Output/protein/tmp/${Candidate_ORFs} > tmp && mv tmp ./Output/protein/tmp/${Candidate_ORFs}
		rm -R "${tBLASTn_candidates}.transdecoder_dir"

		# echo "Creating peptide BLAST database for ${Species} ${Reference_gene}"
		makeblastdb -dbtype prot -in "./Output/protein/tmp/${Candidate_ORFs}" -out "./Output/protein/tmp/${Candidate_ORFs_db}"


		BLASTp_outputfile="${Species}_${Reference_gene}_BLASTp_output.fasta"
		BLASTp_output_IDs="${Species}_${Reference_gene}_BLASTp_IDs.fasta"
		BLASTp_candidates="${Species}_${Reference_gene}_prot.fasta"
		BLASTp_unclear_ORFs="${Species}_${Reference_gene}_unclear_ORFs.fasta"
		BLASTp_missing_ORFs="${Species}_${Reference_gene}_missing_ORFs.fasta"
		BLASTp_all_unclear_ORFs="002_${Reference_gene}_all_unclear_ORFs_IDs.txt"
		BLASTp_all_missing_ORFs="003_${Reference_gene}_all_missing_ORFs_IDs.txt"
		BLASTp_all_missing_ORFs_seqs="003_${Reference_gene}_all_missing_ORFs_nuc_seqs.txt"
		BLASTp_all_candidates="001_${Reference_gene}_all_prot_candidates.fasta"


# ---------- e-value adjustment 2 ----------
	# this is the second place where the minimum e-value for the search should be adjusted 
	# although it is not as important, as long as this one is the same or below the first adjustment
		blastp -query "${gene}" -db "./Output/protein/tmp/${Candidate_ORFs_db}" -task blastp -evalue 1e-1 -out "./Output/protein/${Reference_gene}/BLASTp_output/${BLASTp_outputfile}" -outfmt 6
		
		mkdir "./Output/protein/${Reference_gene}/tmp1"
		cut -f 2 "./Output/protein/${Reference_gene}/BLASTp_output/${BLASTp_outputfile}" > "./Output/protein/${Reference_gene}/tmp1/${BLASTp_outputfile}"
		
		mkdir "./Output/protein/${Reference_gene}/tmp2"
		sort "./Output/protein/${Reference_gene}/tmp1/${BLASTp_outputfile}" > "./Output/protein/${Reference_gene}/tmp2/${BLASTp_outputfile}"
		rm -R "./Output/protein/${Reference_gene}/tmp1"

		uniq "./Output/protein/${Reference_gene}/tmp2/${BLASTp_outputfile}" > "./Output/protein/${Reference_gene}/BLASTp_output/IDs/${BLASTp_output_IDs}"
		rm -R "./Output/protein/${Reference_gene}/tmp2"

		./Apps/fastagrep.pl -f "./Output/protein/${Reference_gene}/BLASTp_output/IDs/${BLASTp_output_IDs}" -X "./Output/protein/tmp/${Candidate_ORFs}" > "./Output/protein/${Reference_gene}/${BLASTp_candidates}"
		
		sed -i "s,_y_y_y_y_p,.p,g" "./Output/protein/${Reference_gene}/BLASTp_output/IDs/${BLASTp_output_IDs}"
		sed -i "s,_y_y_y_y_p,.p,g" "./Output/protein/${Reference_gene}/BLASTp_output/${BLASTp_outputfile}"
		sed -i "s,_y_y_y_y_p,.p,g" "./Output/protein/${Reference_gene}/${BLASTp_candidates}"
		sed -i "s,_y_y_y_y_p,.p,g" "./Output/nucleotide/${Reference_gene}/${tBLASTn_candidates}"

		mkdir "./Output/protein/${Reference_gene}/tmp3"
		mkdir "./Output/protein/${Reference_gene}/tmp4"
		mkdir "./Output/protein/${Reference_gene}/tmp5"
		mkdir "./Output/protein/${Reference_gene}/tmp6"
		mkdir "./Output/protein/${Reference_gene}/tmp7"
		mkdir "./Output/protein/${Reference_gene}/tmp8"

		cp "./Output/protein/${Reference_gene}/BLASTp_output/IDs/${BLASTp_output_IDs}" "./Output/protein/${Reference_gene}/tmp3/${BLASTp_output_IDs}"
		sed -i "s,_x_x_x_x_p,\t,g" "./Output/protein/${Reference_gene}/tmp3/${BLASTp_output_IDs}"
		cut -f 1 "./Output/protein/${Reference_gene}/tmp3/${BLASTp_output_IDs}" > "./Output/protein/${Reference_gene}/tmp4/${BLASTp_output_IDs}"
		sort "./Output/protein/${Reference_gene}/tmp4/${BLASTp_output_IDs}" > "./Output/protein/${Reference_gene}/tmp5/${BLASTp_output_IDs}"
		uniq -d "./Output/protein/${Reference_gene}/tmp5/${BLASTp_output_IDs}" > "./Output/protein/${Reference_gene}/unclear_ORFs/${BLASTp_unclear_ORFs}"


		uniq "./Output/protein/${Reference_gene}/tmp5/${BLASTp_output_IDs}" > "./Output/protein/${Reference_gene}/tmp6/${BLASTp_output_IDs}"
		cat "./Output/nucleotide/${Reference_gene}/tBLASTn_output/IDs/${tBLASTn_output_IDs}" "./Output/protein/${Reference_gene}/tmp6/${BLASTp_output_IDs}" > "./Output/protein/${Reference_gene}/tmp7/${BLASTp_missing_ORFs}"
		sort "./Output/protein/${Reference_gene}/tmp7/${BLASTp_missing_ORFs}" > "./Output/protein/${Reference_gene}/tmp8/${BLASTp_missing_ORFs}"
		uniq -u "./Output/protein/${Reference_gene}/tmp8/${BLASTp_missing_ORFs}" > "./Output/protein/${Reference_gene}/missing_ORFs/${BLASTp_missing_ORFs}"

		rm -R "./Output/protein/${Reference_gene}/tmp3"
		rm -R "./Output/protein/${Reference_gene}/tmp4"
		rm -R "./Output/protein/${Reference_gene}/tmp5"
		rm -R "./Output/protein/${Reference_gene}/tmp6"
		rm -R "./Output/protein/${Reference_gene}/tmp7"
		rm -R "./Output/protein/${Reference_gene}/tmp8"

		sed -i "s,_x_x_x_x_p,.p,g" "./Output/protein/${Reference_gene}/BLASTp_output/IDs/${BLASTp_output_IDs}"
		sed -i "s,_x_x_x_x_p,.p,g" "./Output/protein/${Reference_gene}/BLASTp_output/${BLASTp_outputfile}"
		sed -i "s,_x_x_x_x_p,.p,g" "./Output/protein/${Reference_gene}/${BLASTp_candidates}"


		./Apps/fastagrep.pl -f "./Output/protein/${Reference_gene}/missing_ORFs/${BLASTp_missing_ORFs}" -X "./Output/nucleotide/${Reference_gene}/${tBLASTn_candidates}" > "./Output/protein/${Reference_gene}/missing_ORFs/tmp/${BLASTp_missing_ORFs}"


	done

	# create a list of all missing ORFs that could not be translated
	cat ./Output/protein/${Reference_gene}/missing_ORFs/tmp/*.fasta > "./Output/protein/${Reference_gene}/${BLASTp_all_missing_ORFs_seqs}"
	rm -R "./Output/protein/${Reference_gene}/missing_ORFs/tmp"

	# create a list of all ORFs that where unclear (if a BLASt hit was detected in more than 1 ORF of a single sequence)
	cat ./Output/protein/${Reference_gene}/unclear_ORFs/*.fasta > "./Output/protein/${Reference_gene}/${BLASTp_all_unclear_ORFs}"
	cat ./Output/protein/${Reference_gene}/missing_ORFs/*.fasta > "./Output/protein/${Reference_gene}/${BLASTp_all_missing_ORFs}"

	rm -R "./Output/protein/${Reference_gene}/unclear_ORFs"
	rm -R "./Output/protein/${Reference_gene}/missing_ORFs"
	
	cat ./Output/protein/${Reference_gene}/*.fasta > "./Output/protein/${Reference_gene}/${BLASTp_all_candidates}"


	# use sed command to append the gene name as the first line to the list of all_prot_candidates
	headline="Precursor candidates: ${Reference_gene}"
	output1="tmp1_${Reference_gene}.fasta"
	sed "1s/^/\n\n\n\n${headline} \n\n/" Output/protein/${Reference_gene}/*_all_prot_candidates.fasta > "Output/tmp/${output1}"

done
rm -R "./Output/protein/tmp"
rm Transcriptomes/*.nhr
rm Transcriptomes/*.nin
rm Transcriptomes/*.nsq

# create a single file with all detected genes from all transcriptomes. List is separated by genes.
cat ./Output/tmp/*.fasta > "./Output/protein/All_neuropeptide_prot_candidates.fasta"
rm -R "./Output/tmp"

echo " "
echo "finished"
echo " "
echo "candidate sequences are in the folder 'Output' "
echo "sequences that were detected by tBLASTn but could not get translated are in 'Output/protein/.../003_..._all_missing_ORFs_nuc_seqs.txt'"
