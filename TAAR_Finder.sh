#!/bin/bash


#SBATCH --job-name=TAAR_On   # Job name

eval "$(conda shell.bash hook)"
conda activate olfactory


LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load R/4.2.0-foss-2021a
module load BLAST/2.12.0-Linux_x86_64
module load EMBOSS/6.2.0-goolf-1.7.20
module load SAMtools/1.15-GCC-10.3.0
module load MAFFT/7.467-GCCcore-7.3.0-with-extensions
module load IQ-TREE/2.0-rc1-foss-2018b
module load Python/3.9.5-GCCcore-10.3.0
module load FASTX-Toolkit/0.0.14-goolf-1.7.20


dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"

####Initalize arguments.



genome=$1
TAAR_database=$2
blast_database=$3
scripts_folder_location=$4 ; scripts_location=`echo "$scripts_folder_location" | sed 's/\/$//'`
maximum_intron_length=$5
number_of_thread=$6
multiple_exon_search=$7  ##TRUE OR FALSE. If the organism had no 2 exons TAAR, then the results should be similar, excepted a faster running time for the single exon mode
evalue=$8
tm_filter=$9

maximum_intron_length_half=$((maximum_intron_length / 2))


###Makeblastdb so we can blast genes against the genome

if test -f "$genome.ndb" ; then echo "Genome blast database already exist" ; else makeblastdb -in $genome -dbtype nucl ; fi 
if test -f "$genome.fai" ; then echo "Genome fai file already exist" ; else samtools faidx $genome ; fi 


####Perform a tblastn with known TAAR genes against the genome, with an evalue of 1e-5

tblastn -query $TAAR_database -db $genome -evalue $evalue -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out TAAR_vs_Genome.blastn -num_threads $number_of_thread 


cut -f1 TAAR_vs_Genome.blastn  | sort | uniq -c | sed 's/^ *//g'| sed 's/ /,/g' | awk -F, ' ($1 < 3000) ' | cut -f2 -d "," > good_tblastn_id ; for i in `cat good_tblastn_id` ; do grep "^$i	" TAAR_vs_Genome.blastn >> Filtered_TAAR_vs_Genome.blastn ; done 
mv Filtered_TAAR_vs_Genome.blastn TAAR_vs_Genome.blastn

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Extraction of single exon TAAR genes  #########################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


####Extract non-overlapping best-hit sequences (first round of tblastn with the evalue of 1e-5) and extend 1000bp upstream and 1000bp downstream

Rscript $scripts_location/R_script_numero1.R 


###For each best-hits we will extract the corresponding sequence

xargs samtools faidx $genome < Best_hits_filtered.tsv > Best_hits_filtered.fasta


####Rename sequences to remove the ":" that is not recognized by emboss getorf
sed -i 's/:/-/g' Best_hits_filtered.fasta


####Extract ORF that are atleast 850bp in size and rename output sequences so that the start coordinate is the start of the ORF (same for end)
getorf -sequence Best_hits_filtered.fasta -outseq orf_list.fasta -minsize 850 -find 3



###Rename results of getorf in order to have good fasta headers

sed -i 's/(REVERSE SENSE)/_reverse/g' orf_list.fasta
grep ">" orf_list.fasta | sed 's/>//g' > oldfastaheaders
sed 's/-/	/g' oldfastaheaders | sed 's/_[0-9] \[/	/g' | sed 's/\]//g' | sed 's/_reverse/	reverse/g'  > ren_oldfastaheaders

IFS=$'\n' #treat the file line by line in the following for loop

for line in `cat ren_oldfastaheaders` ; do
		scaffold=`echo "$line" | cut -f1`

		if grep -q "reverse" <<< "$line" ; then 
			coord_start=`echo "$line" | awk 'BEGIN {FS="	"} {sum+=$2+$5-1} END {print sum}'`
			coord_end=`echo "$line" | awk 'BEGIN {FS="	"} {sum+=$2+$4-1} END {print sum}'`
		else
			coord_start=`echo "$line" | awk 'BEGIN {FS="	"} {sum+=$2+$4-1} END {print sum}'`
			coord_end=`echo "$line" | awk 'BEGIN {FS="	"} {sum+=$2+$5-1} END {print sum}'`

		fi

		echo "$scaffold-$coord_start-$coord_end" 


done > newfastaheaders

paste -d "\t" oldfastaheaders newfastaheaders > renaming_file

perl $scripts_location/rename_fasta.pl renaming_file orf_list.fasta > renamed_orf_list.fasta


###Due to the -1000/+1000 extension, some orf could be found twice so we remove identical sequences

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' renamed_orf_list.fasta | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > renamed_orf_list.fasta_uniq.fa

###Lets translate DNA to Prot

transeq renamed_orf_list.fasta_uniq.fa renamed_orf_list.fasta_uniq.prot ; sed -i 's/_1$//g' renamed_orf_list.fasta_uniq.prot


###We apply a first filter to ORFs : blastp sequences against a database contining OR, TAAR, V2R, V1R , T1R, T2R, and other GPCR genes

blastp -query renamed_orf_list.fasta_uniq.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads $number_of_thread

###Retain only genes that best match to TAAR
grep -i "TAAR-Receptor" blastp_result | cut -f1 | sort | uniq > TAAR_from_blast.list
xargs samtools faidx renamed_orf_list.fasta_uniq.prot < TAAR_from_blast.list > Fasta_TAAR_from_blast.prot


###Second gene filter : Lets align our sequences and perform a ML tree with known TAAR genes (max 200 optimization round)
mafft --add Fasta_TAAR_from_blast.prot --keeplength $scripts_location/TAAR_plus_TAARL_plus_outgroup.aln > Putative_taar_plus_known_taar_plus_outgroup.prot.aln
iqtree -s Putative_taar_plus_known_taar_plus_outgroup.prot.aln -st AA -nt $number_of_thread -m JTT+F+I+G4 -redo -n 200

###We use a script to remove outgroup sequences that are not TAAR

grep ">" $scripts_location/Outgroup_seqs.prot | sed 's/>//g' > Outgroups_id.txt
grep ">" $scripts_location/taar-l_cd_hit_70.prot | sed 's/>//g' > taar_l.id
grep ">" $scripts_location/TAAR_database_reformat_cdhit_80.aln | sed 's/>//g' > TAAR_known.id
Rscript $scripts_location/Tree_Parser_TAAR.R #Generate the file "Current_species_TAAR.txt" 

###Remaining sequences are TAAR
xargs samtools faidx renamed_orf_list.fasta_uniq.fa < Current_species_TAAR.txt > Functionnal_TAARs_multifasta_singleexon.fa


###Remove sequences with N (=ambigous sequences)
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_TAARs_multifasta_singleexon.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Functionnal_TAARs_multifasta_singleexon.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_TAARs_multifasta_singleexon.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > unclear_Functionnal_TAARs_multifasta_singleexon.fa


#Translate complete sequences
transeq clear_Functionnal_TAARs_multifasta_singleexon.fa clear_Functionnal_TAARs_multifasta_singleexon.prot ; sed -i 's/_1$//g' clear_Functionnal_TAARs_multifasta_singleexon.prot



###Extract the coordinates of functionnals ORs found. Then tblastn with functional sequences found

grep ">" Functionnal_TAARs_multifasta_singleexon.fa | sed 's/>//g' | sed 's/-/	/g' > Coordinates_Functionnal_TAAR.txt
if [ `wc -l < Coordinates_Functionnal_TAAR.txt` -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_Functionnal_TAAR.txt ; fi



tblastn -query clear_Functionnal_TAARs_multifasta_singleexon.prot -db $genome -evalue $evalue -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out tblastn_functionnal_taar_vs_genome_current_species.tblastn -num_threads $number_of_thread

#Remove problematic sequences
cut -f1 tblastn_functionnal_taar_vs_genome_current_species.tblastn  | sort | uniq -c | sed 's/^ *//g'| sed 's/ /,/g' | awk -F, ' ($1 < 4000) ' | cut -f2 -d "," > good_tblastn_id ; for i in `cat good_tblastn_id` ; do grep "$i" tblastn_functionnal_taar_vs_genome_current_species.tblastn >> Filtered_tblastn_functionnal_taar_vs_genome_current_species.tblastn ; done
cat TAAR_vs_Genome.blastn Filtered_tblastn_functionnal_taar_vs_genome_current_species.tblastn > tblastn_functionnal_and_known_taar_vs_genome.tblastn


cat $TAAR_database clear_Functionnal_TAARs_multifasta_singleexon.prot > Complete_TAAR_db.prot


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Extraction of multiple exon TAAR genes - First round ##########################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


if [ $multiple_exon_search == "TRUE" ] ; then
	### Now we will look for potential multiple exon genes.
	
	
	#First , let's filter regions to identify those which may contain multiple exon genes. 
	#The script remove regions with already identified genes and merge blast hits if there is <=Xbp between them (with X = half the max intron length)
	
	Rscript $scripts_location/R_sript_merge_multiexon.R $maximum_intron_length_half
	
	#extract regions containing putative multiple exon genes
	
	xargs samtools faidx $genome < Potential_multiple_exon_regions.tsv > Potential_multiple_exon_regions.fa
	
	
	#Retain only besthit that have a TAAR gene in atleast the 3 first hits
	blastx -query Potential_multiple_exon_regions.fa -db $blast_database -max_target_seqs 3 -outfmt '6 qseqid sseqid' -out blastx_blast_regions.tsv -num_threads $number_of_thread
	grep "TAAR-Receptor" blastx_blast_regions.tsv | cut -f1 | sort | uniq > TAAR_best_hits_regions.tsv
	xargs samtools faidx $genome < TAAR_best_hits_regions.tsv > TAAR_best_hits_regions.fa
	
	
	#Launch exonerate with identified species TAAR and with the TAAR database against regions that may contain multiple exon TAAR genes
	
	#split the database in X to run exonerate on X threads
	mkdir Splitted_db
	$scripts_location/exonerate-2.2.0-x86_64/bin/fastasplit -f Complete_TAAR_db.prot -c 30 --output Splitted_db
	
	mkdir Exonerate_raw_results_folder
	
	for i in Splitted_db/* ; do
		file_name=`echo $i | sed 's/Splitted_db\///g'`
		sbatch -W -c 4 --qos=6hours --wrap="$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --ryo '%tcs' --minintron 50 --maxintron $maximum_intron_length Splitted_db/$file_name TAAR_best_hits_regions.fa > Exonerate_raw_results_folder/$file_name.exo.rslt ; sleep 10" &
	done
	
	echo "Exonerate running -- Wait"
	
	wait
	
	echo "Exonerate research done"
	
	
	#Concatenate all exonerate results in one file
	
	cat Exonerate_raw_results_folder/*.exo.rslt > Exonerate_results.txt
	
	
	
	#extract vulgar lines
	
	grep "vulgar" Exonerate_results.txt > vulgar_lines.txt 
	
	
	#Extract interesting columns of vulgar lines
	
	sed 's/vulgar: //g' vulgar_lines.txt | cut -f1,2,3,5,6,7,8,9 -d " " > vulgar_lines_parsed.txt
	
	
	#Count the number of introns using vulgar lines
	
	awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/ I /,"") "\t" NR}' vulgar_lines.txt > number_introns_per_line.txt
	grep -v "count" number_introns_per_line.txt | cut -f1  > intron_numbers.txt 
	
	#Add the intron number to vulgar lines 
	
	paste -d " " vulgar_lines_parsed.txt intron_numbers.txt > vulgar_lines_intron_numbers.txt
	
	
	#Add informations about the best blastp results of each exonerate predicted genes
	
	
	sed -n '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/p' Exonerate_results.txt | sed 's/C4 Alignment:.*//g' | sed 's/Hostname:.*//g' | sed 's/Command line:.*//g' | sed 's/^--$//g' | sed 's/-- completed exonerate analysis.*//g' | sed 's/# --- END OF GFF DUMP ---//g' | sed 's/^#$/>seq_to_rename/g' > List_exonerate_cds.fasta #extract all predicted genes sequences
	transeq List_exonerate_cds.fasta List_exonerate_cds.prot #translate sequences
	sed 's/ /_/g' vulgar_lines_intron_numbers.txt > sequences_names.txt #extract exonerate vulgar line to rename sequences
	awk '/^>/ { printf("%s_%s\n",$0,i++);next;} { print $0;}' List_exonerate_cds.prot > List_exonerate_cds_renamed.prot #first round of rename
	grep ">" List_exonerate_cds_renamed.prot | sed 's/>//g' > old_names.txt #extract names
	paste -d "\t" old_names.txt sequences_names.txt > renaming_file #crate a file for rename_fasta.pl
	perl $scripts_location/rename_fasta.pl renaming_file List_exonerate_cds_renamed.prot > List_exonerate_cds.prot #completely rename sequences with the exonerate vulgar line
	
	#Perform a blastp against all known TAAR genes 
	blastp -query List_exonerate_cds.prot -db $scripts_location/TAAR_plus_TAARL_database.prot -outfmt '6 qseqid sseqid evalue' -out all_blastp.txt -max_target_seqs 1 -num_threads $number_of_thread
	
	#Extract the best blastp hit for each predicted cds
	[ -e all_blastp_parsed.txt ] && rm all_blastp_parsed.txt
	for i in `cat sequences_names.txt` ; do if grep -q "$i" all_blastp.txt ; then grep -m1 "$i" all_blastp.txt | cut -f2,3 >> all_blastp_parsed.txt ; else echo "NoQuery	99999" >> all_blastp_parsed.txt ; fi ; done 
	sed -i 's/	/ /g' all_blastp_parsed.txt
	paste -d " " vulgar_lines_intron_numbers.txt all_blastp_parsed.txt > vulgar_lines_intron_numbers_blastrslt.txt
	
	
	
	#Now let's filter exonerate results. The Rscript below retain non-overlapping best-hits from the exonerate file. I applies several filters : Number of exon >= number of blast hit in the same region + Filter on query length followed by filter on scaffold length
	
	
	Rscript $scripts_location/R_script_exonerate_parsing_new.R #Generate the file "Filtered_Potential_multiple_exon_regions.tsv"
	
	
	#Transform the coordinates to match scaffold coordinates and add the query 
	
	
	
	IFS=$'\n'
	
	for line in `cat Filtered_Potential_multiple_exon_regions.tsv` ; do
		
		query=`echo "$line" | cut -f7`
		scaffold=`echo "$line" | cut -f1`
		hit_start=`echo "$line" | cut -f2`
		hit_end=`echo "$line" | cut -f3`
	
	
		echo "$scaffold:$hit_start-$hit_end	$query"
	
	done > Correct_coordinates_for_exonerate.tsv
	
	
	
	#Let's predict genes on our regions 
	
	IFS=$'\n'
	
	
	mkdir Multiple_exon_predictions
	
	
	for line in `cat Correct_coordinates_for_exonerate.tsv` ; do 
		
		scaffold_s_e=`echo "$line" | cut -f1`
		query=`echo "$line" | cut -f2`
		scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'`
	
		samtools faidx $genome $scaffold_s_e > scaffold.fa
		sed -i 's/:/-/g' scaffold.fa
		samtools faidx $scripts_location/TAAR_plus_TAARL_database.prot $query > query.prot
	
		$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Multiple_exon_predictions/$scaffold_s_e_n.exonerate
		cp scaffold.fa Multiple_exon_predictions/$scaffold_s_e_n.fasta
	
	done
	
	
	
	### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
	
	
	#Result folder
	mkdir Filtered_predictions
	
	
	for file in Multiple_exon_predictions/*.exonerate ; do
	
		#extract some infos from file name
		file_name=`echo "$file" | sed 's/.*\///g'`
		file_name_reduced=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
		fasta_file_name=`echo "$file_name" | sed 's/exonerate/fasta/g'`
		initial_header=`grep ">" Multiple_exon_predictions/$fasta_file_name | sed 's/>//g'`
	
	
	
		#Test if the predicted gene is an OR or not
		awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
		transeq predicted_cds.fa predicted_cds.prot
		blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 10
		
	
		#Lets continue only if the best match is an OR
		#if grep -q -i "olfactory\|odorant" blastp_result ; then 
		if grep -q -i "TAAR-Receptor" blastp_result ; then
	
			#Define the scaffold  name
			scaffold=`echo "$file" | sed 's/.*\///g' | sed 's/-.*//g'`
			#Define the strand on which the predicted gene is
			strand=`grep "	similarity	" $file | cut -f7`
			#Define the first position of the query on the target sequence
			first_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
			#Define the last position of the query on the target sequence
			second_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
			
			
			#Lets extract CDS if the gene is on the negative strand
			if [ $strand == "-" ] ; then 
			
			
				#If strand is minus, then the first position is:
				target_end=$((first_hit_range + 1))
				#And we will went to extend this by 500bp to be sure to have the potentiel start codon
				target_extanded_end=$((first_hit_range + 500))
			
				#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
				samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
				samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
			
				#remove fasta header of extanded region files
				grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
				grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
			
				#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
	
				grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
				for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
				grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
	
				#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
				cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
				sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
				getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 850 -find 3
				if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then name_seq=`grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/ .*//g' | sed 's/>//g'` ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $name_seq > mygoodseq.orf ; mv mygoodseq.orf Filtered_predictions/$file_name_reduced.ORF  ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
			
	
	
				#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
				if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
			
					transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
					$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Multiple_exon_predictions/$fasta_file_name > verif_coord.exo
					if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Multiple_exon_predictions/$fasta_file_name > verif_coord.exo ; fi
	
					extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
					cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
					cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
					exon_number=`grep "	exon	" verif_coord.exo | wc -l`
					sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
			
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					#number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | wc -l`
					number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
					if [ "$exon_number" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
	
				#If not ORF found, then determinate the gene state
			
				elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
			
					stop_codon_state="FALSE"
					edge_state="FALSE"
					frameshift_state="FALSE"
			
					##Stop codon checking
			
					#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
					transeq predicted_cds.fa predicted_cds.prot
			
					#Estimate the interval on which we wil search stop codons. 
					query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
					query_total_length=`grep -m1 "$query_name" $scripts_location/TAAR_plus_TAARL_database.prot.fai | cut -f2`
					query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
			
					#Lets see if we find stop codon before the five_percent_position
					stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
			
					if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
			
			
					##Frameshift checking
			
					#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
					#We remove border exons if there are less than 60nt in length. Run as iteration.
			
					grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
					awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
			
					#Check for the presence of frameshift
					frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
					if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
					if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
			
			
			
					##Edge checking
			
					#Check if the gene is at a conting border
					#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
					gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
					gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
			
					extracted_scaffold_start=`grep ">" Multiple_exon_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
			
					true_start_coord=$((extracted_scaffold_start + gene_start_coord))
					true_end_coord=$((extracted_scaffold_start + gene_end_coord))
			
					#First check if these coordinates are near the end of scaffolds (<100 bp)
					if [ "$true_start_coord" -le '100' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
					scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
					diff_lengths=$((scaffold_length - true_end_coord))
					if [ "$diff_lengths" -le '100' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
					
					#Now check if there are consecutive N near the gene that could indicate conting end
					extanded_start_coord=$((true_start_coord - 200))
					extanded_end_coord=$((true_end_coord + 200))
			
					#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
					consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
					if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
					if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
			
			
					##Extract the sequence
			
					[ -e Current_exon_rev.txt ] && rm Current_exon_rev.txt
			
					#Extract the corresponding sequence
					for line in `cat Correct_exons.txt` ; do
						start_pos=`echo "$line" | cut -f1`
						end_pos=`echo "$line" | cut -f2`
			
						samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
						revseq Current_exon.fa Current_exon_rev.fa
			
						
			
						#add the reversed sequence to a text file
						grep -v ">" Current_exon_rev.fa >> Current_exon_rev.txt
			
					done
			
					#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
					exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
			
					header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
					sed -e "1i>$header_name\\" Current_exon_rev.txt > Filtered_predictions/$file_name_reduced.PSEU
					sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
					cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
					sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					#number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | wc -l`
					number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
					if [ "$exon_nb" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
					if [ "$exon_nb" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
	
				fi
	
	
	
			#Lets make the same steps with slight modifications for the + strand
			elif [ $strand == "+" ] ; then 
			
				#If strand is minus, then the first position is:
				target_end=$((second_hit_range + 1))
				#And we will went to extend this by 500bp to be sure to have the potentiel start codon
				target_extanded_end=$((second_hit_range + 500))
			
				#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
				samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
				samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
			
				#remove fasta header of extanded region files
				grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
				grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
			
				#Extract the CDS sequence predicted by exonerate and remove fasta header
	
				grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
				for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
				grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
			
				#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
				cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
				sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
				getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 850 -find 3 -reverse FALSE
			
			
				#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
				if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
			
					transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
					$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Multiple_exon_predictions/$fasta_file_name > verif_coord.exo
					if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Multiple_exon_predictions/$fasta_file_name > verif_coord.exo ; fi
	
					extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
					cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
					cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
					exon_number=`grep "	exon	" verif_coord.exo | wc -l`
					sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
				
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					#number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | wc -l`
					number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
					if [ "$exon_number" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
						
				#If not ORF found, then determinate the gene state
			
				elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
			
					stop_codon_state="FALSE"
					edge_state="FALSE"
					frameshift_state="FALSE"
			
				##Stop codon checking
			
					#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
					transeq predicted_cds.fa predicted_cds.prot
			
					#Estimate the interval on which we wil search stop codons. 
					query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
					query_total_length=`grep -m1 "$query_name" $scripts_location/TAAR_plus_TAARL_database.prot.fai | cut -f2`
					query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
			
					#Lets see if we find stop codon before the five_percent_position
					stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
			
					if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
			
			
					##Frameshift checking
			
					#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
					#We remove border exons if there are less than 60nt in length. Run as iteration.
			
					grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
					awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
			
					#Check for the presence of frameshift
					frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
					if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
					if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
			
			
			
					##Edge checking
			
					#Check if the gene is at a conting border
					#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
					gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
					gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
			
					extracted_scaffold_start=`grep ">" Multiple_exon_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
			
					true_start_coord=$((extracted_scaffold_start + gene_start_coord))
					true_end_coord=$((extracted_scaffold_start + gene_end_coord))
			
					#First check if these coordinates are near the end of scaffolds (<100 bp)
					if [ "$true_start_coord" -le '100' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
					scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
					diff_lengths=$((scaffold_length - true_end_coord))
					if [ "$diff_lengths" -le '100' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
					
					#Now check if there are consecutive N near the gene that could indicate conting end
					extanded_start_coord=$((true_start_coord - 200))
					extanded_end_coord=$((true_end_coord + 200))
			
					#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
					consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
					if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
					if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
			
			
					##Extract the sequence
			
					[ -e Current_exon.txt ] && rm Current_exon.txt
					#Extract the corresponding sequence
					for line in `cat Correct_exons.txt` ; do
						start_pos=`echo "$line" | cut -f1`
						end_pos=`echo "$line" | cut -f2`
			
						samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
			
						#add the reversed sequence to a text file
						grep -v ">" Current_exon.fa >> Current_exon.txt
			
					done
			
			
					#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
					exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
			
					header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
					sed -e "1i>$header_name\\" Current_exon.txt > Filtered_predictions/$file_name_reduced.PSEU
					sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
					cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
					sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					#number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | wc -l`
					number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi
					if [ "$exon_nb" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
					if [ "$exon_nb" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
	
	
				fi
			fi
		fi
	
	done
	
	
	
	
	
	
	
	
	
	#Now that we have filtered all our results, we can concatenate the results
	
	for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
	cat Filtered_predictions/*.ORF > Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
	cat Filtered_predictions/*.PSEU > Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
	cat Filtered_predictions/*.CDSP > Frameshift_less_Pseudogenes.fa
	
	
	
	
	
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################## Extraction of multiple exon TAAR genes - Second round  ########################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	
	
	rm -r Multiple_exon_predictions/
	rm -r Filtered_predictions/
	
	
	#Extract coordinates of found gene in the first round 
	
	grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_MulitpleExon_TAAR.txt
	grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_MulitpleExon_TAAR.txt
	
	
	if [ `wc -l < Coordinates_MulitpleExon_TAAR.txt` -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_MulitpleExon_TAAR.txt ; fi
		
	#Now let's filter a second time exonerate results. The Rscript below retain non-overlapping best-hits from the exonerate file. 
	#I applies several filters : Number of exon >= number of blast hit in the same region + Filter on query length followed by filter on scaffold length
	
	
	Rscript $scripts_location/R_script_exonerate_parsing_new_second.R #Generate the file "Filtered_Potential_multiple_exon_regions.tsv"
	
	
	#Transform the coordinates to match scaffold coordinates and add the query 
	
	
	
	IFS=$'\n'
	
	for line in `cat Filtered_Potential_multiple_exon_regions.tsv` ; do
		
		query=`echo "$line" | cut -f7`
		scaffold=`echo "$line" | cut -f1`
		hit_start=`echo "$line" | cut -f2`
		hit_end=`echo "$line" | cut -f3`
	
	
		echo "$scaffold:$hit_start-$hit_end	$query"
	
	done > Correct_coordinates_for_exonerate.tsv
	
	
	
	#Let's predict genes on our regions 
	
	IFS=$'\n'
	
	
	mkdir Multiple_exon_predictions
	
	
	for line in `cat Correct_coordinates_for_exonerate.tsv` ; do 
		
		scaffold_s_e=`echo "$line" | cut -f1`
		query=`echo "$line" | cut -f2`
		scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'`
	
		samtools faidx $genome $scaffold_s_e > scaffold.fa
		sed -i 's/:/-/g' scaffold.fa
		samtools faidx $scripts_location/TAAR_plus_TAARL_database.prot $query > query.prot
	
		$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Multiple_exon_predictions/$scaffold_s_e_n.exonerate
		cp scaffold.fa Multiple_exon_predictions/$scaffold_s_e_n.fasta
	
	done
	
	
	
	### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
	
	
	#Result folder
	mkdir Filtered_predictions
	
	
	for file in Multiple_exon_predictions/*.exonerate ; do
	
		#extract some infos from file name
		file_name=`echo "$file" | sed 's/.*\///g'`
		file_name_reduced=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
		fasta_file_name=`echo "$file_name" | sed 's/exonerate/fasta/g'`
		initial_header=`grep ">" Multiple_exon_predictions/$fasta_file_name | sed 's/>//g'`
	
	
	
		#Test if the predicted gene is an OR or not
		awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
		transeq predicted_cds.fa predicted_cds.prot
		blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 10
		
	
		#Lets continue only if the best match is an OR
		#if grep -q -i "olfactory\|odorant" blastp_result ; then 
		if grep -q -i "TAAR-Receptor" blastp_result ; then
	
			#Define the scaffold  name
			scaffold=`echo "$file" | sed 's/.*\///g' | sed 's/-.*//g'`
			#Define the strand on which the predicted gene is
			strand=`grep "	similarity	" $file | cut -f7`
			#Define the first position of the query on the target sequence
			first_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
			#Define the last position of the query on the target sequence
			second_hit_range=`grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
			
			
			#Lets extract CDS if the gene is on the negative strand
			if [ $strand == "-" ] ; then 
			
			
				#If strand is minus, then the first position is:
				target_end=$((first_hit_range + 1))
				#And we will went to extend this by 500bp to be sure to have the potentiel start codon
				target_extanded_end=$((first_hit_range + 500))
			
				#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
				samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
				samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
			
				#remove fasta header of extanded region files
				grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
				grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
			
				#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
	
				grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
				for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
				grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
	
				#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
				cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
				sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
				getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 850 -find 3
				if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then name_seq=`grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/ .*//g' | sed 's/>//g'` ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $name_seq > mygoodseq.orf ; mv mygoodseq.orf Filtered_predictions/$file_name_reduced.ORF  ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
			
	
	
				#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
				if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
			
					transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
					$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Multiple_exon_predictions/$fasta_file_name > verif_coord.exo
					if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Multiple_exon_predictions/$fasta_file_name > verif_coord.exo ; fi
	
					extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
					cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
					cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
					exon_number=`grep "	exon	" verif_coord.exo | wc -l`
					sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					#number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | wc -l`
					number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
					if [ "$exon_number" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
	
				#If not ORF found, then determinate the gene state
			
				elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
			
					stop_codon_state="FALSE"
					edge_state="FALSE"
					frameshift_state="FALSE"
			
					##Stop codon checking
			
					#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
					transeq predicted_cds.fa predicted_cds.prot
			
					#Estimate the interval on which we wil search stop codons. 
					query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
					query_total_length=`grep -m1 "$query_name" $scripts_location/TAAR_plus_TAARL_database.prot.fai | cut -f2`
					query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
			
					#Lets see if we find stop codon before the five_percent_position
					stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
			
					if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
			
			
					##Frameshift checking
			
					#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
					#We remove border exons if there are less than 60nt in length. Run as iteration.
			
					grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
					awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
			
					#Check for the presence of frameshift
					frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
					if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
					if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
			
			
			
					##Edge checking
			
					#Check if the gene is at a conting border
					#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
					gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
					gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
			
					extracted_scaffold_start=`grep ">" Multiple_exon_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
			
					true_start_coord=$((extracted_scaffold_start + gene_start_coord))
					true_end_coord=$((extracted_scaffold_start + gene_end_coord))
			
					#First check if these coordinates are near the end of scaffolds (<100 bp)
					if [ "$true_start_coord" -le '100' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
					scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
					diff_lengths=$((scaffold_length - true_end_coord))
					if [ "$diff_lengths" -le '100' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
					
					#Now check if there are consecutive N near the gene that could indicate conting end
					extanded_start_coord=$((true_start_coord - 200))
					extanded_end_coord=$((true_end_coord + 200))
			
					#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
					consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
					if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
					if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
			
			
					##Extract the sequence
			
					[ -e Current_exon_rev.txt ] && rm Current_exon_rev.txt
			
					#Extract the corresponding sequence
					for line in `cat Correct_exons.txt` ; do
						start_pos=`echo "$line" | cut -f1`
						end_pos=`echo "$line" | cut -f2`
			
						samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
						revseq Current_exon.fa Current_exon_rev.fa
			
						
			
						#add the reversed sequence to a text file
						grep -v ">" Current_exon_rev.fa >> Current_exon_rev.txt
			
					done
			
					#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
					exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
			
					header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
					sed -e "1i>$header_name\\" Current_exon_rev.txt > Filtered_predictions/$file_name_reduced.PSEU
					sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
					cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
					sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					#number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | wc -l`
					number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi  
					if [ "$exon_nb" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
					if [ "$exon_nb" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
	
	
	
				fi
	
	
	
			#Lets make the same steps with slight modifications for the + strand
			elif [ $strand == "+" ] ; then 
			
				#If strand is minus, then the first position is:
				target_end=$((second_hit_range + 1))
				#And we will went to extend this by 500bp to be sure to have the potentiel start codon
				target_extanded_end=$((second_hit_range + 500))
			
				#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
				samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
				samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
			
				#remove fasta header of extanded region files
				grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
				grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
			
				#Extract the CDS sequence predicted by exonerate and remove fasta header
	
				grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
				for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
				grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
			
				#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
				cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
				sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
				getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 850 -find 3 -reverse FALSE
			
			
				#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
				if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
			
					transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
					$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Multiple_exon_predictions/$fasta_file_name > verif_coord.exo
					if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Multiple_exon_predictions/$fasta_file_name > verif_coord.exo ; fi
	
	
					extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
					cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
					cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
					exon_number=`grep "	exon	" verif_coord.exo | wc -l`
					sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
				
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					#number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | wc -l`
					number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
					if [ "$exon_number" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
						
				#If not ORF found, then determinate the gene state
			
				elif [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -lt 1 ] ; then 
			
					stop_codon_state="FALSE"
					edge_state="FALSE"
					frameshift_state="FALSE"
			
				##Stop codon checking
			
					#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
					transeq predicted_cds.fa predicted_cds.prot
			
					#Estimate the interval on which we wil search stop codons. 
					query_name=`grep "Query: " $file | sed 's/.*Query: //g'`
					query_total_length=`grep -m1 "$query_name" $scripts_location/TAAR_plus_TAARL_database.prot.fai | cut -f2`
					query_start_position=`grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
			
					#Lets see if we find stop codon before the five_percent_position
					stop_codon_nb=`grep -v ">" predicted_cds.prot | fold -c1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l` #number of stop codons before the ten percent pos
			
					if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
			
			
					##Frameshift checking
			
					#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
					#We remove border exons if there are less than 60nt in length. Run as iteration.
			
					grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
					awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
			
					#Check for the presence of frameshift
					frameshift_nb=`grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}'`
					if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
					if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
			
			
			
					##Edge checking
			
					#Check if the gene is at a conting border
					#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
					gene_start_coord=`cut -f1 Correct_exons.txt | sort -n | head -1`
					gene_end_coord=`cut -f2 Correct_exons.txt | sort -n | tail -1`
			
					extracted_scaffold_start=`grep ">" Multiple_exon_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
			
					true_start_coord=$((extracted_scaffold_start + gene_start_coord))
					true_end_coord=$((extracted_scaffold_start + gene_end_coord))
			
					#First check if these coordinates are near the end of scaffolds (<100 bp)
					if [ "$true_start_coord" -le '100' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
					scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
					diff_lengths=$((scaffold_length - true_end_coord))
					if [ "$diff_lengths" -le '100' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
					
					#Now check if there are consecutive N near the gene that could indicate conting end
					extanded_start_coord=$((true_start_coord - 200))
					extanded_end_coord=$((true_end_coord + 200))
			
					#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
					consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
					if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
					if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
			
			
					##Extract the sequence
			
					[ -e Current_exon.txt ] && rm Current_exon.txt
					#Extract the corresponding sequence
					for line in `cat Correct_exons.txt` ; do
						start_pos=`echo "$line" | cut -f1`
						end_pos=`echo "$line" | cut -f2`
			
						samtools faidx Multiple_exon_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
			
						#add the reversed sequence to a text file
						grep -v ">" Current_exon.fa >> Current_exon.txt
			
					done
			
			
					#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
					exon_nb=`wc -l Correct_exons.txt | sed 's/ .*//g'`
			
					header_name=`echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
					sed -e "1i>$header_name\\" Current_exon.txt > Filtered_predictions/$file_name_reduced.PSEU
					sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
					cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
					sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
					#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
					samtools faidx $genome $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
					makeblastdb -in Verification_scaffold.fa -dbtype nucl
					#number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | wc -l`
					number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi
					if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
					if [ "$exon_nb" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
					if [ "$exon_nb" -gt 2 ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
	
	
	
				fi
			fi
		fi
	
	done
	
	
	
	
	
	
	
	
	#Now that we have filtered all our results, we can concatenate the results
	
	for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
	cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
	cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
	cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa
	
	
	
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################## Verification of mulitple exon genes  #######################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	
	
	
	#Remove sequences with identical headers
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_multiple_exon.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Pseudogenes_multiple_exon_uniq.fa ; mv Pseudogenes_multiple_exon_uniq.fa Pseudogenes_multiple_exon.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Frameshift_less_Pseudogenes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Frameshift_less_Pseudogenes_uniq.fa ; mv Frameshift_less_Pseudogenes_uniq.fa Frameshift_less_Pseudogenes.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Potential_multiple_exon_CDS.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Potential_multiple_exon_CDS_uniq.fa ; mv Potential_multiple_exon_CDS_uniq.fa Potential_multiple_exon_CDS.fa
	
	
	#Verify complete genes with a blast
	transeq Potential_multiple_exon_CDS.fa Potential_multiple_exon_CDS.prot ; sed -i 's/_1$//g' Potential_multiple_exon_CDS.prot 
	blastp -query Potential_multiple_exon_CDS.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out second_blastp_result -max_target_seqs 1 -num_threads $number_of_thread
	grep -i "TAAR-Receptor" second_blastp_result | cut -f1 | sort | uniq > multiple_exon_TAAR_from_blast.list
	xargs samtools faidx Potential_multiple_exon_CDS.prot < multiple_exon_TAAR_from_blast.list > Filtered_Potential_multiple_exon_CDS.prot 
	
	
	#verification of pseudogenes with a blastp (with the file with frameshifts removed)
	blastp -query Frameshift_less_Pseudogenes.fa -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out pseudo_blastp_result -max_target_seqs 1 -num_threads $number_of_thread
	grep -i "TAAR-Receptor" pseudo_blastp_result | cut -f1 | sort | uniq > multiple_exon_pseudos_from_blast.list
	fasta_formatter -i Frameshift_less_Pseudogenes.fa > test.fasta ; mv test.fasta Frameshift_less_Pseudogenes.fa
	fasta_formatter -i Pseudogenes_multiple_exon.fa > test.fasta ; mv test.fasta Pseudogenes_multiple_exon.fa
	xargs samtools faidx Frameshift_less_Pseudogenes.fa < multiple_exon_pseudos_from_blast.list > Verified_Pseudogenes_multiple_exon_frameshift_less.fa
	
	#Remove ambigous sequences
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Verified_Pseudogenes_multiple_exon_frameshift_less.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\X/)' | tr "\t" "\n" > Verified_Pseudogenes_multiple_exon_frameshift_less_clean.prot


	## Lets make a tree
	nb_seq=`grep -c ">" Verified_Pseudogenes_multiple_exon_frameshift_less.fa`
	if [ "$nb_seq" -gt "0" ] ; then
		mafft --add Verified_Pseudogenes_multiple_exon_frameshift_less.fa --keeplength $scripts_location/TAAR_plus_TAARL_plus_outgroup.aln > Multiple_exon_verification.aln
	else
		cp $scripts_location/TAAR_plus_TAARL_plus_outgroup.aln ./Multiple_exon_verification.aln
	fi
	
	#also add functionnal genes and make the tree
	
	nb_seq=`grep -c ">" Filtered_Potential_multiple_exon_CDS.prot`
	if [ "$nb_seq" -gt "0" ] ; then
		mafft --add Filtered_Potential_multiple_exon_CDS.prot --keeplength Multiple_exon_verification.aln > Putative_taar_plus_known_taar_plus_outgroup.prot.aln
		iqtree -s Putative_taar_plus_known_taar_plus_outgroup.prot.aln -st AA -nt $number_of_thread -m JTT+F+I+G4 -redo -n 200
	else
		cp Multiple_exon_verification.aln Putative_taar_plus_known_taar_plus_outgroup.prot.aln
		iqtree -s Putative_taar_plus_known_taar_plus_outgroup.prot.aln -st AA -nt $number_of_thread -m JTT+F+I+G4 -redo -n 200
	fi
	
	
	#Parse the tree and retain good sequences
	rm Current_species_TAAR.txt
	Rscript $scripts_location/Tree_Parser_TAAR.R 
	grep "FALSE\|TRUE" Current_species_TAAR.txt  > TAAR_pseudos.id
	grep -v "FALSE\|TRUE" Current_species_TAAR.txt  > TAAR_functionnal.id
	xargs samtools faidx Potential_multiple_exon_CDS.fa < TAAR_functionnal.id > Verified_Potential_multiple_exon_CDS.fa
	xargs samtools faidx Pseudogenes_multiple_exon.fa < TAAR_pseudos.id > Verified_Pseudogenes_TAAR_genes.fa 
	
	
	#remove pseudogenes with a length < 200nt (150 for single exon pseudogene later, but for two exons, put 200)
	seqkit seq -m 200 Verified_Pseudogenes_TAAR_genes.fa  > Verified_Pseudogenes_TAAR_genes_s200.fa ; mv Verified_Pseudogenes_TAAR_genes_s200.fa Verified_Pseudogenes_TAAR_genes.fa
	
	
	#Add the coordinates of found genes to the file "Coordinates_Functionnal_TAAR.txt"
	
	grep ">" Verified_Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_Functionnal_TAAR.txt
	grep ">" Verified_Pseudogenes_TAAR_genes.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_Functionnal_TAAR.txt
	
	
	#remove ambigous sequences in complete genes 
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Verified_Potential_multiple_exon_CDS.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Verified_Potential_multiple_exon_CDS.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Verified_Potential_multiple_exon_CDS.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > unclear_Verified_Potential_multiple_exon_CDS.fa
	

	#No need to clean remaining pseudogenes as this file was already filtered for ambigous sequences before making the tree
	cp Verified_Pseudogenes_TAAR_genes.fa clear_Verified_Pseudogenes_TAAR_genes.fa


fi

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Extraction of singe exon TAAR pseudogenes  ####################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

#Launch a Rscript. This will find regions corresponding to single exon 
#pseudogenes/ truncated genes by parsing the tblastn results and removing already analyzed regions

Rscript $scripts_location/R_script_numero2_2.R #This Rscript generate "Pseudo_truncated_coordinates.tsv" with 4 columns : scaffold, coord_start, coord_end, query, non-extanded start, non-extanded stop

##Coordinates_Functionnal_TAAR ??
#extract regions containing putative TAAR pseudogenes

cut -f1 Pseudo_truncated_coordinates.tsv > scaffolds.txt
cut -f2,3 Pseudo_truncated_coordinates.tsv | sed 's/	/-/g' > coordinates.txt
paste -d ":" scaffolds.txt coordinates.txt > Pseudogenes_regions.tsv

xargs samtools faidx $genome < Pseudogenes_regions.tsv > Pseudogenes_regions.fa


#Retain only besthit that have a TAAR gene in atleast the 3 first hits
blastx -query Pseudogenes_regions.fa -db $blast_database -max_target_seqs 3 -outfmt '6 qseqid sseqid' -out blastx_blast_regions.tsv -num_threads $number_of_thread
grep "TAAR-Receptor" blastx_blast_regions.tsv | cut -f1 | sort | uniq > TAAR_best_hits_regions.tsv
xargs samtools faidx $genome < TAAR_best_hits_regions.tsv > TAAR_best_hits_regions.fa


IFS=$'\n'
for line in `cat TAAR_best_hits_regions.tsv` ; do 
	scaffold=`echo "$line" | sed 's/:/	/g' | sed 's/-/	/g'| cut -f1`
	start=`echo "$line" | sed 's/:/	/g' | sed 's/-/	/g'| cut -f2`
	stop=`echo "$line" | sed 's/:/	/g' | sed 's/-/	/g'| cut -f3`

	grep "$scaffold.*$start.*$stop" Pseudo_truncated_coordinates.tsv >> Pseudo_truncated_coordinates_filtered.tsv
done


## remove regions less than 150bp ##

IFS=$'\n'
for line in `cat Pseudo_truncated_coordinates_filtered.tsv` ; do
	scaffold=`echo $line | cut -f1`
	true_start=`echo $line | cut -f5` 
	true_end=`echo $line | cut -f6` 

	length=$((true_end - true_start))

	if [ "$length" -ge '150' ] ; then echo "$line" >> Filtered_regions.tsv ; fi 

done

mv Filtered_regions.tsv Pseudo_truncated_coordinates_filtered.tsv



#Lets loop over filtered regions to characterize genes

IFS=$'\n'

[ -e Pseudogenes_single_exon.fa ] && rm Pseudogenes_single_exon.fa
[ -e Frameshiftless_Pseudogenes_single_exon.prot ] && rm Frameshiftless_Pseudogenes_single_exon.prot



for line in `cat Pseudo_truncated_coordinates_filtered.tsv` ; do

	#initialize values
	scaffold=`echo $line | cut -f1`
	start=`echo $line | cut -f2` #extanded start
	stop=`echo $line | cut -f3` #extanded stop
	query=`echo $line | cut -f4`
	scaffold_length=`grep "^$scaffold	" $genome.fai | cut -f2`
	true_start=`echo $line | cut -f5` #real best-hit start coord
	true_end=`echo $line | cut -f6` #real best-hit end coord
	strand=`echo $line | cut -f7`



	#initialize states
	stop_codon_state="FALSE"
	edge_state="FALSE"
	frameshift_state="FALSE"


	#Extract the genomic region as well as the best query that matched on it with blast
	samtools faidx $genome $scaffold:$start-$stop > Genomic_region.fa
	makeblastdb -in Genomic_region.fa -dbtype nucl  #make a blast database with the genomic region 
	samtools faidx Complete_TAAR_db.prot $query > Gene.prot

	blastx -query Genomic_region.fa -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastx_result -max_target_seqs 1 -num_threads 10
	 
	if grep -q -i "TAAR-Receptor" blastx_result ; then

		#tblast the query against the region and print the target sequence in the blast result (without gaps to detect frameshits)
		tblastn -query Gene.prot -db Genomic_region.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe sseq' > tblast_result_stop.tsv
		tblastn -gapopen 32767 -gapextend 32767 -query Gene.prot -db Genomic_region.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe sseq' > tblast_result_fs.tsv
	
		#The column 12 correspond to the different detected frames
		cut -f12 tblast_result_fs.tsv > Frame_detection.txt
	
		#The column 13 correspond to the target sequence, usefull to detect stop codon (asterix)
		cut -f13 tblast_result_stop.tsv > Stop_detection.txt
	
	
		#Check the number of hsp in the tblast result
		number_hsp=`cat Frame_detection.txt | wc -l` 
	
		number_diff_frames=`cat Frame_detection.txt | sort | uniq | wc -l`
	
	
		#if [ "$number_diff_frames" -ge '1' ] ; then frameshift_state="TRUE" ; fi
		if [ "$number_diff_frames" -gt '1' ] ; then frameshift_state="TRUE" ; fi
	
	
		#Check the number of stop codons present in the HSPs
		number_stops=`grep -o "\*" Stop_detection.txt | wc -l`
	
		if [ "$number_stops" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
	
	
		#Check if the gene is at a conting border (as we did for multiple exon gene search)
	
		if [ "$true_start" -le '100' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
		diff_lengths=$((scaffold_length - true_end))
		if [ "$diff_lengths" -le '100' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
		
		extanded_start_coord=$((start - 100))
		extanded_end_coord=$((stop + 100))
		consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
		if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
		if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
	
	
		#Extract sequences and rename fasta header to get the full infos
	
	
		if [ $strand == "+" ] ; then
	
			samtools faidx $genome $scaffold:$true_start-$true_end > temporary_rslt.fa
			grep -v ">" temporary_rslt.fa > temporary_rslt.txt
			header_name=`echo "$scaffold-$true_start-$true_end---1_exon-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
			sed -e "1i>$header_name\\"  temporary_rslt.txt > temporary_rslt_renamed.fa
			sed '/^[[:space:]]*$/d' temporary_rslt_renamed.fa >> Pseudogenes_single_exon.fa
	
			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 Gene.prot Genomic_region.fa > Exonerate_single_exon_pseudo
			awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' Exonerate_single_exon_pseudo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
			transeq predicted_cds.fa predicted_cds.prot ; rm predicted_cds.fa
			sed "s/>.*/>$header_name/g" predicted_cds.prot >> Frameshiftless_Pseudogenes_single_exon.prot
	
	
	
	
	
		else
	
			samtools faidx $genome $scaffold:$true_start-$true_end > temporary_rslt.fa
			revseq temporary_rslt.fa temporary_rslt_rev.fa
			grep -v ">" temporary_rslt_rev.fa > temporary_rslt_rev.txt
			header_name=`echo "$scaffold-$true_start-$true_end---1_exon-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'`
			sed -e "1i>$header_name\\" temporary_rslt_rev.txt > temporary_rslt_renamed.fa
			sed '/^[[:space:]]*$/d' temporary_rslt_renamed.fa >> Pseudogenes_single_exon.fa
	
			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 Gene.prot Genomic_region.fa > Exonerate_single_exon_pseudo 
			awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' Exonerate_single_exon_pseudo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
			transeq predicted_cds.fa predicted_cds.prot ; rm predicted_cds.fa
			sed "s/>.*/>$header_name/g" predicted_cds.prot >> Frameshiftless_Pseudogenes_single_exon.prot
	
	
		fi
	fi
done	




#IFS=$'\n' ; for line in `cat Pseudo_truncated_coordinates.tsv` ; do scaffold=`echo $line | cut -f1` ; start=`echo $line | cut -f2` ; stop=`echo $line | cut -f3` ; query=`echo $line | cut -f4` ; scaffold_length=`grep "^$scaffold " $genome.fai | cut -f2` ; true_start=`echo $line | cut -f5` ; true_end=`echo $line | cut -f6` ; strand=`echo $line | cut -f7` ; stop_codon_state="FALSE" ; edge_state="FALSE" ; frameshift_state="FALSE" ; samtools faidx $genome $scaffold:$start-$stop > Genomic_region.fa ; makeblastdb -in Genomic_region.fa -dbtype nucl ; samtools faidx Complete_TAAR_db.prot $query > Gene.prot ; blastx -query Genomic_region.fa -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastx_result -max_target_seqs 1 -num_threads 10 ; if grep -q -i "TAAR-Receptor" blastx_result ; then tblastn -query Gene.prot -db Genomic_region.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe sseq' > tblast_result_stop.tsv ; tblastn -gapopen 32767 -gapextend 32767 -query Gene.prot -db Genomic_region.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe sseq' > tblast_result_fs.tsv ; cut -f12 tblast_result_fs.tsv > Frame_detection.txt ; cut -f13 tblast_result_stop.tsv > Stop_detection.txt ; number_hsp=`cat Frame_detection.txt | wc -l`  ; number_diff_frames=`cat Frame_detection.txt | sort | uniq | wc -l` ; if [ "$number_diff_frames" -gt '1' ] ; then frameshift_state="TRUE" ; fi ; number_stops=`grep -o "\*" Stop_detection.txt | wc -l` ; if [ "$number_stops" -ge '1' ] ; then stop_codon_state="TRUE" ; fi ; if [ "$true_start" -le '100' ] ; then edge_state="TRUE" ; fi  ; diff_lengths=$((scaffold_length - true_end)) ; if [ "$diff_lengths" -le '100' ] ; then edge_state="TRUE" ; fi ; extanded_start_coord=$((start - 100)) ; extanded_end_coord=$((stop + 100)) ; consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1` ; if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi ; if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi  ; if [ $strand == "+" ] ; then samtools faidx $genome $scaffold:$true_start-$true_end > temporary_rslt.fa ; grep -v ">" temporary_rslt.fa > temporary_rslt.txt ; header_name=`echo "$scaffold-$true_start-$true_end---1_exon-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'` ; sed -e "1i>$header_name\\"  temporary_rslt.txt > temporary_rslt_renamed.fa ; sed '/^[[:space:]]*$/d' temporary_rslt_renamed.fa >> Pseudogenes_single_exon.fa ; $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 Gene.prot Genomic_region.fa > Exonerate_single_exon_pseudo ; awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' Exonerate_single_exon_pseudo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa ; transeq predicted_cds.fa predicted_cds.prot ; rm predicted_cds.fa ; sed "s/>.*/>$header_name/g" predicted_cds.prot >> Frameshiftless_Pseudogenes_single_exon.prot ; else samtools faidx $genome $scaffold:$true_start-$true_end > temporary_rslt.fa ; revseq temporary_rslt.fa temporary_rslt_rev.fa ; grep -v ">" temporary_rslt_rev.fa > temporary_rslt_rev.txt ; header_name=`echo "$scaffold-$true_start-$true_end---1_exon-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'` ; sed -e "1i>$header_name\\" temporary_rslt_rev.txt > temporary_rslt_renamed.fa ; sed '/^[[:space:]]*$/d' temporary_rslt_renamed.fa >> Pseudogenes_single_exon.fa ; $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 Gene.prot Genomic_region.fa > Exonerate_single_exon_pseudo  ; awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' Exonerate_single_exon_pseudo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa ; transeq predicted_cds.fa predicted_cds.prot ; rm predicted_cds.fa ; sed "s/>.*/>$header_name/g" predicted_cds.prot >> Frameshiftless_Pseudogenes_single_exon.prot ; fi ;fi ; done

# Remove ambigous sequences 

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Frameshiftless_Pseudogenes_single_exon.prot | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\X/)' | tr "\t" "\n" > Frameshiftless_Pseudogenes_single_exon_clean.prot
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Frameshiftless_Pseudogenes_single_exon.prot | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\X/)' | tr "\t" "\n" > Frameshiftless_Pseudogenes_single_exon_unclean.prot
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_single_exon.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > Pseudogenes_single_exon_clean.fa 
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_single_exon.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > Pseudogenes_single_exon_unclean.fa 


# Filter sequences with a length <150nt

seqkit seq -m 150 Pseudogenes_single_exon_clean.fa  > Pseudogenes_single_exon_clean_s200.fa 
fasta_formatter -i Pseudogenes_single_exon_clean_s200.fa > test.fasta ; mv test.fasta Pseudogenes_single_exon_clean_s200.fa
grep ">" Pseudogenes_single_exon_clean_s200.fa | sed 's/>//g' > retained_id ; xargs samtools faidx Frameshiftless_Pseudogenes_single_exon.prot < retained_id > Frameshiftless_Pseudogenes_single_exon_s200.prot


#Verify pseudogenes with a ML tree

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_single_exon_clean_s200.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Pseudogenes_single_exon_uniq.fa ; mv Pseudogenes_single_exon_uniq.fa Pseudogenes_single_exon_clean_s200.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Frameshiftless_Pseudogenes_single_exon_s200.prot | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Frameshiftless_Pseudogenes_single_exon_uniq.fa ; mv Frameshiftless_Pseudogenes_single_exon_uniq.fa Frameshiftless_Pseudogenes_single_exon_s200.prot


nb_seq=`grep -c ">" Frameshiftless_Pseudogenes_single_exon_s200.prot`
if [ "$nb_seq" -gt "0" ] ; then
	mafft --add Frameshiftless_Pseudogenes_single_exon_s200.prot --keeplength $scripts_location/TAAR_plus_TAARL_plus_outgroup.aln > Putative_taar_plus_known_taar_plus_outgroup.prot.aln
	iqtree -s Putative_taar_plus_known_taar_plus_outgroup.prot.aln -st AA -nt $number_of_thread -m JTT+F+I+G4 -redo -n 200
	Rscript $scripts_location/Tree_Parser_TAAR.R 

	xargs samtools faidx Pseudogenes_single_exon_clean_s200.fa < Current_species_TAAR.txt > Verified_Pseudogenes_single_exon_clean_s200.fa

else

	echo "" > Verified_Pseudogenes_single_exon_clean_s200.fa


fi



##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

#Add the number of exon to functional genes
grep ">" clear_Functionnal_TAARs_multifasta_singleexon.fa | sed 's/>//g' > single_exon_id.txt
sed -e 's/$/---1_exons/' single_exon_id.txt > single_exon_id_edit.txt
paste -d "\t" single_exon_id.txt single_exon_id_edit.txt > renaming_single_exon.txt
perl $scripts_location/rename_fasta.pl renaming_single_exon.txt clear_Functionnal_TAARs_multifasta_singleexon.fa > temporary.fasta ; mv temporary.fasta clear_Functionnal_TAARs_multifasta_singleexon.fa


#Merge pseudogenes with multiple and one exon

cat clear_Verified_Pseudogenes_TAAR_genes.fa Verified_Pseudogenes_single_exon_clean_s200.fa > Final_Pseudogenes.fa

#Merge complete genes with multiple and one exon
cat clear_Verified_Potential_multiple_exon_CDS.fa clear_Functionnal_TAARs_multifasta_singleexon.fa  > Combined_Functionnal_TAAR.fa


#Check if complete genes have 7tm domain determined with phobius or TMHMM


#First check with phobius
transeq Combined_Functionnal_TAAR.fa Combined_Functionnal_TAAR.prot ; sed -i 's/_1$//g' Combined_Functionnal_TAAR.prot #translate CDS
perl $scripts_location/phobius/phobius.pl -long Combined_Functionnal_TAAR.prot > Phobius_verification.txt #run phonius in long mode
grep ">" Combined_Functionnal_TAAR.prot | sed 's/>//g' > gene_id.txt #extract cds id

for gene in `cat gene_id.txt` ; do nb_transm=`sed '/'"$gene"'/,/\/\//!d;/\/\//q' Phobius_verification.txt | grep "TRANSMEM" | wc -l` ; echo "$gene,$nb_transm" ; done > Gene_NbTm.tsv
awk 'BEGIN{FS=",";OFS=","}($2>=7){print $1;}' Gene_NbTm.tsv > Phobius_genes_with_7tm.txt
awk 'BEGIN{FS=",";OFS=","}($2<7){print $1;}' Gene_NbTm.tsv > Phobius_genes_without_7tm.txt

#Now, with TMHMM

$scripts_location/tmhmm-2.0c/bin/tmhmm Combined_Functionnal_TAAR.prot > tmhmm_verification.txt
grep "Number of predicted TMHs:" tmhmm_verification.txt | sed 's/# //g' | sed 's/ Number of predicted TMHs:  /,/g' | awk -F "," '{ if(($2 >= 7)) { print $1} }' > tmhmm_genes_with_7tm.txt
grep "Number of predicted TMHs:" tmhmm_verification.txt | sed 's/# //g' | sed 's/ Number of predicted TMHs:  /,/g' | awk -F "," '{ if(($2 < 7)) { print $1} }' > tmhmm_genes_without_7tm.txt 


#Combine results of predictions
cat Phobius_genes_with_7tm.txt tmhmm_genes_with_7tm.txt | sort | uniq > Genes_with_7tm.txt
sort gene_id.txt > gene_id_sorted.txt ; sort Genes_with_7tm.txt > sorted_Genes_with_7tm.txt
comm -3 gene_id_sorted.txt sorted_Genes_with_7tm.txt > Genes_without_7tm.txt


#Print the result in the fasta file

for gene in `cat gene_id.txt` ; do
	pred_phobius="FALSE"
	pred_tmhmm="FALSE"

	if grep -q "$gene" Phobius_genes_with_7tm.txt ; then pred_phobius="TRUE" ; fi
	if grep -q "$gene" tmhmm_genes_with_7tm.txt ; then pred_tmhmm="TRUE" ; fi 

	if [ $pred_phobius == "TRUE" ] && [ $pred_tmhmm == "TRUE" ] ; then 
		new_gene_name="$gene---phobius-tmhmm"
	elif [ $pred_phobius == "TRUE" ] && [ $pred_tmhmm == "FALSE" ] ; then 
		new_gene_name="$gene---phobius" 
	elif [ $pred_phobius == "FALSE" ] && [ $pred_tmhmm == "TRUE" ] ; then
		new_gene_name="$gene---tmhmm"
	else
		new_gene_name="$gene"
	fi 

	echo "$new_gene_name" 

done > New_gene_name_with_predictions.txt


paste -d "\t" gene_id.txt New_gene_name_with_predictions.txt > renaming_file_tm.txt

perl $scripts_location/rename_fasta.pl renaming_file_tm.txt Combined_Functionnal_TAAR.fa > temporary.fasta ; mv temporary.fasta Combined_Functionnal_TAAR.fa




#Merge all files that contain ambigous sequences

cat Pseudogenes_single_exon_unclean.fa unclear_Functionnal_TAARs_multifasta_singleexon.fa unclear_Verified_Potential_multiple_exon_CDS.fa > Ambigous_TAAR.fasta


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Correction of 1 exon predictions  ###########################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

if [ $multiple_exon_search == "TRUE" ] ; then
	## For TAAR genes, it is not unusual that 1 exon gene found in the first steps of the script are in fact two exons genes that we did not find complete
	## Most of the time, these genes will not have a predicted GPCR domain by tmhmm or phobius
	## We will verify that and correct bad predictions
	
	#First lets check 1 exon genes without predicted 7tm
	
	#grep id of such genes and take their coordinates
	grep ">" Combined_Functionnal_TAAR.fa | grep "1_exons" | sed 's/>//g' | grep -v "phobius" | grep -v "tmhmm" | sort > id_one_exon_verification.txt
	sed 's/-/	/g' id_one_exon_verification.txt | cut -f1,2,3 > Coordinates_genes_verifications.tsv
	
	
	if [ `wc -l < Coordinates_genes_verifications.tsv` -ge 1 ] ; then 
	
		#grep id of all good genes and their coordinates as well as pseudogenes
		grep ">" Combined_Functionnal_TAAR.fa | sed 's/>//g' | sort > rest
		comm -23 rest id_one_exon_verification.txt > Good_genes.txt
		sed 's/-/	/g' Good_genes.txt | cut -f1,2,3 > Coordinates_correct_genes.tsv
		
		grep ">" Final_Pseudogenes.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_correct_genes.tsv
		grep ">" Ambigous_TAAR.fasta | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_correct_genes.tsv
		
		
		#Now lets launch a Rscript to extend the regions of genes (max 3000bp upstream and downstream) to be verificated but do not extend over good genes or found pseudogenes
		
		
		Rscript $scripts_location/Extend_predicted_one_exon.R $maximum_intron_length
		
		
		xargs samtools faidx $genome < Extended_regions.tsv > Extended_regions.fa
		
		
		mkdir Exonerate_second_row_raw_results_folder
		
		
		for i in Splitted_db/* ; do
			file_name=`echo $i | sed 's/Splitted_db\///g'`
			sbatch -W -c 4 --qos=6hours --wrap="$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Splitted_db/$file_name Extended_regions.fa > Exonerate_second_row_raw_results_folder/$file_name.exo.rslt ; sleep 10" &
		done
		
		echo "Exonerate running -- Wait"
		
		wait
		
		echo "Exonerate research done"
		
		
		#Concatenate all exonerate results in one file
		
		cat Exonerate_second_row_raw_results_folder/*.exo.rslt > Exonerate_second_results.txt
		
		
		#extract vulgar lines
		
		grep "vulgar" Exonerate_second_results.txt > vulgar_lines.txt 
		
		
		#extract interesting columns of vulgar lines
		#query, query_start, query_end, scaffold, scaffold_start, scaffold_end, strand, exonerate_score
		
		sed 's/vulgar: //g' vulgar_lines.txt | cut -f1,2,3,5,6,7,8,9 -d " " > vulgar_lines_parsed.txt
		
		
		#count the number of introns using vulgar lines
		
		
		IFS=$'\n'
		
		awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/ I /,"") "\t" NR}' vulgar_lines.txt > number_introns_per_line.txt
		grep -v "count" number_introns_per_line.txt | cut -f1  > intron_numbers.txt 
		
		
		#add the intron number to vulgar lines 
		
		paste -d " " vulgar_lines_parsed.txt intron_numbers.txt > vulgar_lines_intron_numbers.txt
		
		
		sed -n '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/p' Exonerate_second_results.txt | sed 's/C4 Alignment:.*//g' | sed 's/Hostname:.*//g' | sed 's/Command line:.*//g' | sed 's/^--$//g' | sed 's/-- completed exonerate analysis.*//g' | sed 's/# --- END OF GFF DUMP ---//g' | sed 's/^#$/>seq_to_rename/g' > List_exonerate_cds.fasta #extract all predicted genes sequences
		transeq List_exonerate_cds.fasta List_exonerate_cds.prot #translate sequences
		sed 's/ /_/g' vulgar_lines_intron_numbers.txt > sequences_names.txt #extract exonerate vulgar line to rename sequences
		awk '/^>/ { printf("%s_%s\n",$0,i++);next;} { print $0;}' List_exonerate_cds.prot > List_exonerate_cds_renamed.prot #first round of rename
		grep ">" List_exonerate_cds_renamed.prot | sed 's/>//g' > old_names.txt #extract names
		paste -d "\t" old_names.txt sequences_names.txt > renaming_file #crate a file for rename_fasta.pl
		perl $scripts_location/rename_fasta.pl renaming_file List_exonerate_cds_renamed.prot > List_exonerate_cds.prot #completely rename sequences with the exonerate vulgar line
		
		#Perform the blastp
		blastp -query List_exonerate_cds.prot -db $scripts_location/TAAR_plus_TAARL_database.prot -outfmt '6 qseqid sseqid evalue' -out all_blastp.txt -max_target_seqs 1 -num_threads $number_of_thread
		
		#Extract the information
		[ -e all_blastp_parsed.txt ] && rm all_blastp_parsed.txt
		for i in `cat sequences_names.txt` ; do if grep -q "$i" all_blastp.txt ; then grep -m1 "$i" all_blastp.txt | cut -f2,3 >> all_blastp_parsed.txt ; else echo "NoQuery	99999" >> all_blastp_parsed.txt ; fi ; done 
		sed -i 's/	/ /g' all_blastp_parsed.txt
		paste -d " " vulgar_lines_intron_numbers.txt all_blastp_parsed.txt > vulgar_lines_intron_numbers_blastrslt.txt
		
		
		#Now let's filter exonerate results. The Rscript below retain non-overlapping best-hits from the exonerate file. If two hits are overlapping, then it keep 
		#genes that are >230aa (almost full length to be an OR) and  we keep the one which take the least amount of place on the scaffold. This generate the file "Filtered_Potential_multiple_exon_regions.tsv"
		
		
		Rscript $scripts_location/R_script_exonerate_parsing_new.R
		
		
		IFS=$'\n'
		
		for line in `cat Filtered_Potential_multiple_exon_regions.tsv` ; do
		
			query=`echo "$line" | cut -f7`
			scaffold=`echo "$line" | cut -f1`
			scaff_start=`echo "$line" | cut -f2`
			scaff_end=`echo "$line" | cut -f3`
		
			echo "$scaffold:$scaff_start-$scaff_end	$query"
		
		done > Correct_coordinates_for_exonerate.tsv
		
		
		mkdir Verification_folder
		
		
		
		IFS=$'\n'
		
		[ -e Corrected_predictions.fa ] && rm Corrected_predictions.fa
		
		
		for line in `cat Correct_coordinates_for_exonerate.tsv` ; do
		
			scaffold=`echo "$line" | sed 's/:/	/g' | sed 's/-/	/g' | cut -f1`
			extend_start=`echo "$line" | sed 's/:/	/g' | sed 's/-/	/g' | cut -f2`
			extend_stop=`echo "$line" | sed 's/:/	/g' | sed 's/-/	/g' | cut -f3`
			best_query=`echo "$line" | cut -f2`
		
		
			#Extract the extended regions
			region_name=`echo "$scaffold-$extend_start-$extend_stop"`
			samtools faidx $genome $scaffold:$extend_start-$extend_stop > Verification_folder/$region_name.fasta
			sed -i 's/:/-/g' Verification_folder/$region_name.fasta
			initial_header=`grep ">" Verification_folder/$region_name.fasta | sed 's/>//g'`
		
		
		
			#Extract the known TAAR genes that have the best exonerate score on our region and perform an exhaustive search with this gene
			samtools faidx $scripts_location/TAAR_plus_TAARL_database.prot $best_query > Current_query.prot
			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E TRUE --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 Current_query.prot Verification_folder/$region_name.fasta > Verification_folder/$region_name.exo
		
			#Continue the loop only if predicted gene has more than 1 exon and if the best blastp result is a TAAR gene
			nb_exon=`grep -o "	exon	" Verification_folder/$region_name.exo | wc -l`
			awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' Verification_folder/$region_name.exo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
			transeq predicted_cds.fa predicted_cds.prot
			blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 10
		
		
			if grep -q -i "TAAR-Receptor" blastp_result && [ "$nb_exon" -ge '2' ] ; then
		
		
				#Define the strand
				strand=`grep "	similarity	" Verification_folder/$region_name.exo | cut -f7`
				#Define the first position of the query on the target sequence
				first_hit_range=`grep -m1 "Target range:" Verification_folder/$region_name.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
				#Define the last position of the query on the target sequence
				second_hit_range=`grep -m1 "Target range:" Verification_folder/$region_name.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
			
			
				
				#Lets extract CDS if the gene is on the negative strand
				if [ $strand == "-" ] ; then 
				
					#If strand is minus, then the first position is:
					target_end=$((first_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((first_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Verification_folder/$region_name.fasta $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
					samtools faidx Verification_folder/$region_name.fasta $initial_header:1-$second_hit_range > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
				
					#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
					grep "	exon	"  Verification_folder/$region_name.exo | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					rm Correct_cds.fa ; for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Verification_folder/$region_name.fasta $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Verification_folder/$region_name.ORF -minsize 850 -find 3
					if grep -q -i "reverse" Verification_folder/$region_name.ORF ; then sequence_to_grep=`grep -i "reverse" Verification_folder/$region_name.ORF | sed 's/>//g' | sed 's/ .*//g'`  ; samtools faidx Verification_folder/$region_name.ORF $sequence_to_grep > temporary ; mv temporary Verification_folder/$region_name.ORF ; rm Verification_folder/$region_name.ORF.fai ; else rm Verification_folder/$region_name.ORF ; echo "bad strand" > Verification_folder/$region_name.ORF ; fi
				
					#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
					if [ `grep -c ">" Verification_folder/$region_name.ORF` -ge 1 ] ; then
				
						transeq Verification_folder/$region_name.ORF Verification_folder/$region_name.ORFP
						$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Verification_folder/$region_name.ORFP Verification_folder/$region_name.fasta > verif_coord.exo
						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Verification_folder/$region_name.ORFP Verification_folder/$region_name.fasta > verif_coord.exo ; fi
						extracted_scaffold_start=$extend_start
						cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
						cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=`grep "	exon	" verif_coord.exo | wc -l`
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Verification_folder/$region_name.ORF
				
						#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
						samtools faidx $genome $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
						makeblastdb -in Verification_scaffold.fa -dbtype nucl
						#number_blast_hit=`tblastn -query Verification_folder/$region_name.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | wc -l`
						number_blast_hit=`tblastn -query Verification_folder/$region_name.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`
						if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Verification_folder/$region_name.ORF ; else cat Verification_folder/$region_name.ORF >> Corrected_predictions.fa ; fi 
		
		
		
		
					else
		
						echo "No TAAR gene found"
		
					fi
		
				#Lets make the same steps with slight modifications for the + strand
				elif [ $strand == "+" ] ; then 
				
					#If strand is minus, then the first position is:
					target_end=$((second_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((second_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Verification_folder/$region_name.fasta $initial_header:1-$first_hit_range > Extend_three_prime.fa
					samtools faidx Verification_folder/$region_name.fasta $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
				
					#Extract the CDS sequence predicted by exonerate and remove fasta header
					grep "	exon	"  Verification_folder/$region_name.exo | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Verification_folder/$region_name.fasta $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Verification_folder/$region_name.ORF -minsize 850 -find 3 -reverse FALSE
				
				
					#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
					if [ `grep -c ">" Verification_folder/$region_name.ORF` -ge 1 ] ; then
				
						transeq Verification_folder/$region_name.ORF Verification_folder/$region_name.ORFP
						$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Verification_folder/$region_name.ORFP Verification_folder/$region_name.fasta > verif_coord.exo
						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Verification_folder/$region_name.ORFP Verification_folder/$region_name.fasta > verif_coord.exo ; fi
						extracted_scaffold_start=$extend_start
						cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
						cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=`grep "	exon	" verif_coord.exo | wc -l`
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Verification_folder/$region_name.ORF
					
						#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
						samtools faidx $genome $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
						makeblastdb -in Verification_scaffold.fa -dbtype nucl
						#number_blast_hit=`tblastn -query Verification_folder/$region_name.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | wc -l`
						number_blast_hit=`tblastn -query Verification_folder/$region_name.ORFP -db Verification_scaffold.fa -evalue 1e-20 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`
						if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Verification_folder/$region_name.ORF ; else cat Verification_folder/$region_name.ORF >> Corrected_predictions.fa ; fi 
		
		
		
					else
		
						echo "No TAAR gene found"
		
		
					fi
				fi
			
			else 
		
				echo "No TAAR gene found"
		
		
			fi
		
		done
		
		
		
		
		#lets generate .fai files of potentially bad 1 exon genes and newly predicted genes
		xargs samtools faidx Combined_Functionnal_TAAR.fa < id_one_exon_verification.txt > To_verify_genes_Combined_Functionnal_TAAR.fa
		samtools faidx To_verify_genes_Combined_Functionnal_TAAR.fa
		
		awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Corrected_predictions.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Corrected_predictions.fa 
		awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' clear_Corrected_predictions.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Corrected_predictions.fa
		
		samtools faidx Corrected_predictions.fa
		
		#lets merge these file and use Rscript to keep the best gene (keep 2 exon) when there is overlap
		cat To_verify_genes_Combined_Functionnal_TAAR.fa.fai Corrected_predictions.fa.fai > Combined_old_new_genes.fai
		cut -f1 Combined_old_new_genes.fai | sed 's/---/-/g' | sed 's/-/	/g' | sed 's/_exons//g' > Combined_old_new_genes.tsv
		cut -f1 Combined_old_new_genes.fai > gene_names.txt
		paste -d "\t" Combined_old_new_genes.tsv gene_names.txt > Combined_old_new_genes_to_parse.tsv
		
		
		Rscript $scripts_location/Keep_best_gene.R #Rscript to keep the best gene , generate the file "Genes_to_keep.tsv"
		
		
		cat To_verify_genes_Combined_Functionnal_TAAR.fa Corrected_predictions.fa > ALL_predictions.fa
		xargs samtools faidx ALL_predictions.fa < Genes_to_keep.tsv > Retained_genes.fa
		xargs samtools faidx Combined_Functionnal_TAAR.fa < Good_genes.txt > Good_genes_Combined_Functionnal_TAAR.fa
		
		
		#relaunch phobius and tmhmm of retained genes
		
		
		transeq Retained_genes.fa Retained_genes.prot ; sed -i 's/_1$//g' Retained_genes.prot #translate CDS
		perl $scripts_location/phobius/phobius.pl -long Retained_genes.prot > Phobius_verification.txt #run phonius in long mode
		grep ">" Retained_genes.prot | sed 's/>//g' > gene_id.txt #extract cds id
		
		for gene in `cat gene_id.txt` ; do nb_transm=`sed '/'"$gene"'/,/\/\//!d;/\/\//q' Phobius_verification.txt | grep "TRANSMEM" | wc -l` ; echo "$gene,$nb_transm" ; done > Gene_NbTm.tsv
		awk 'BEGIN{FS=",";OFS=","}($2>=7){print $1;}' Gene_NbTm.tsv > Phobius_genes_with_7tm.txt
		awk 'BEGIN{FS=",";OFS=","}($2<7){print $1;}' Gene_NbTm.tsv > Phobius_genes_without_7tm.txt
		
		#Now, with TMHMM
		
		$scripts_location/tmhmm-2.0c/bin/tmhmm Retained_genes.prot > tmhmm_verification.txt
		grep "Number of predicted TMHs:" tmhmm_verification.txt | sed 's/# //g' | sed 's/ Number of predicted TMHs:  /,/g' | awk -F "," '{ if(($2 >= 7)) { print $1} }' > tmhmm_genes_with_7tm.txt
		grep "Number of predicted TMHs:" tmhmm_verification.txt | sed 's/# //g' | sed 's/ Number of predicted TMHs:  /,/g' | awk -F "," '{ if(($2 < 7)) { print $1} }' > tmhmm_genes_without_7tm.txt 
		
		
		#Print the result in the fasta file
		
		for gene in `cat gene_id.txt` ; do
			pred_phobius="FALSE"
			pred_tmhmm="FALSE"
		
			if grep -q "$gene" Phobius_genes_with_7tm.txt ; then pred_phobius="TRUE" ; fi
			if grep -q "$gene" tmhmm_genes_with_7tm.txt ; then pred_tmhmm="TRUE" ; fi 
		
			if [ $pred_phobius == "TRUE" ] && [ $pred_tmhmm == "TRUE" ] ; then 
				new_gene_name="$gene---phobius-tmhmm"
			elif [ $pred_phobius == "TRUE" ] && [ $pred_tmhmm == "FALSE" ] ; then 
				new_gene_name="$gene---phobius" 
			elif [ $pred_phobius == "FALSE" ] && [ $pred_tmhmm == "TRUE" ] ; then
				new_gene_name="$gene---tmhmm"
			else
				new_gene_name="$gene"
			fi 
		
			echo "$new_gene_name" 
		
		done > New_gene_name_with_predictions.txt
		
		
		paste -d "\t" gene_id.txt New_gene_name_with_predictions.txt > renaming_file_tm.txt
		
		perl $scripts_location/rename_fasta.pl renaming_file_tm.txt Retained_genes.fa > temporary.fasta ; mv temporary.fasta Retained_genes.fa
		
		
		#Merge all functionnal genes
		
		cat Good_genes_Combined_Functionnal_TAAR.fa Retained_genes.fa > Combined_Functionnal_TAAR.fa
	
	fi

fi


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Parsing result files  ######################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################



# If two genes/pseudogenes are overlapping due to the multiple extensions, then keep only the longest found gene/pseudogene

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Ambigous_TAAR.fasta | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Ambigous_TAAR_uniq.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Final_Pseudogenes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Final_Pseudogenes_uniq.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Combined_Functionnal_TAAR.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Combined_Functionnal_TAAR_uniq.fa



nb_seq=`grep -c ">" Ambigous_TAAR_uniq.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" Ambigous_TAAR_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_ambigous_final.tsv
fi

nb_seq=`grep -c ">" Final_Pseudogenes_uniq.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" Final_Pseudogenes_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_pseudogenes_final.tsv
fi

nb_seq=`grep -c ">" Combined_Functionnal_TAAR_uniq.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" Combined_Functionnal_TAAR_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_genes_final.tsv
fi



Rscript $scripts_location/Remove_redundancy.R



IFS=$'\n'


for line in `cat best_genes_functionnal.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Combined_Functionnal_TAAR_uniq.fa | sed 's/>//g' >> functionnal_to_keep.txt ; done
for line in `cat best_genes_ambigous.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Ambigous_TAAR_uniq.fa | sed 's/>//g' >> ambigous_to_keep.txt ; done
for line in `cat best_genes_pseudogenes.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Final_Pseudogenes_uniq.fa | sed 's/>//g' >> pseudogenes_to_keep.txt ; done


xargs samtools faidx Combined_Functionnal_TAAR_uniq.fa < functionnal_to_keep.txt > FINAL_Functionnal_TAAR.fa
xargs samtools faidx Ambigous_TAAR_uniq.fa < ambigous_to_keep.txt > FINAL_Ambigous_TAAR.fa
xargs samtools faidx Final_Pseudogenes_uniq.fa < pseudogenes_to_keep.txt > FINAL_Pseudogenes_TAAR.fa


#Also classify genes without 7tm as pseudogenes


if [ $tm_filter == "TRUE" ] ; then

	grep ">" FINAL_Functionnal_TAAR.fa | grep "phobius\|tmhmm" | sed 's/>//g' > 7tm_genes
	grep ">" FINAL_Functionnal_TAAR.fa | grep -v "phobius\|tmhmm" | sed 's/>//g' > non_7tm_genes
	xargs samtools faidx FINAL_Functionnal_TAAR.fa < 7tm_genes > FINAL_Functionnal_TAAR_7tm.fa 
	xargs samtools faidx FINAL_Functionnal_TAAR.fa < non_7tm_genes >> FINAL_Pseudogenes_TAAR.fa

else
	cp FINAL_Functionnal_TAAR.fa FINAL_Functionnal_TAAR_7tm.fa
fi



## FINAL FILTER

transeq FINAL_Functionnal_TAAR_7tm.fa FINAL_Functionnal_TAAR_7tm.prot ; sed -i 's/_1$//g' FINAL_Functionnal_TAAR_7tm.prot
blastp -query FINAL_Functionnal_TAAR_7tm.prot -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastp_result_custom -max_target_seqs 1 -num_threads $number_of_thread
grep "TAAR" blastp_result_custom | cut -f1 | sort | uniq  > good_sequences
xargs samtools faidx FINAL_Functionnal_TAAR_7tm.fa < good_sequences > temp.fa
mv temp.fa FINAL_Functionnal_TAAR_7tm.fa ; rm good_sequences ; rm *.fai


blastx -query FINAL_Pseudogenes_TAAR.fa -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastx_result_custom -max_target_seqs 1 -num_threads $number_of_thread
grep "TAAR" blastx_result_custom | cut -f1 | sort | uniq  > good_sequences_p
xargs samtools faidx FINAL_Pseudogenes_TAAR.fa < good_sequences_p > temp_p.fa
mv temp_p.fa FINAL_Pseudogenes_TAAR.fa ; rm good_sequences_p ; rm *.fai

blastx -query FINAL_Ambigous_TAAR.fa -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastx_result_custom -max_target_seqs 1 -num_threads $number_of_thread
grep "TAAR" blastx_result_custom | cut -f1 | sort | uniq  > good_sequences_p
xargs samtools faidx FINAL_Ambigous_TAAR.fa < good_sequences_p > temp_p.fa
mv temp_p.fa FINAL_Ambigous_TAAR.fa ; rm good_sequences_p ; rm *.fai



#We have three final file :


nb_functionnal=`grep -c ">" FINAL_Functionnal_TAAR_7tm.fa`
nb_pseudo_edge=`grep -c ">" FINAL_Pseudogenes_TAAR.fa`
nb_ambigous=`grep -c ">" FINAL_Ambigous_TAAR.fa`

echo "Search of TAAR is finished. There are $nb_functionnal potentially functionnal genes, $nb_pseudo_edge pseudogenes or fragments and $nb_ambigous ambigous sequences"

echo "$nb_functionnal	$nb_pseudo_edge	$nb_ambigous" > Results_NbF_NbP_NbA_summary.txt



dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"


#grep -c ">" FINAL_Functionnal_TAAR.fa ; grep -c ">" FINAL_Pseudogenes_TAAR.fa ; grep -c ">" FINAL_Ambigous_TAAR.fa

#If the user want, we can compute a tree with all olfactory receptor found.
#transeq FINAL_Functionnal_TAAR.fa FINAL_Functionnal_TAAR.prot ; sed -i 's/_1$//g' FINAL_Functionnal_TAAR.prot 
#mafft --add FINAL_Functionnal_TAAR.prot --keeplength $scripts_location/TAAR_plus_TAARL_plus_outgroup.aln > RESULTS_Functionnal_genes_plus_outgroup.prot.aln
#iqtree -s RESULTS_Functionnal_genes_plus_outgroup.prot.aln -st AA -nt $number_of_thread -bb 1000 -m JTT+F+I+G4 -redo 

##dups : FINAL_COORD_TAAR_non_retained_similatiry.txt

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## END  #######################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

