#!/bin/bash


#SBATCH --job-name=V1R_Finder   # Job name


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
V1R_database=$2
blast_database=$3
scripts_folder_location=$4 ; scripts_location=`echo "$scripts_folder_location" | sed 's/\/$//'`
maximum_intron_length=$5
number_of_thread=$6
multiple_exon_search=$7
evalue=$8
tm_filter=$9

### Unlike TAAR genes, here if we put TRUE or FALSE in the multiple exon gene search, it will run very different pipelines.
### As TAAR genes, the FALSE pipeline run way faster...
### TRUE is adapted to actinopterygi
### FALSE is adapted to tetrapods


###Makeblastdb so we can blast genes against the genome


if test -f "$genome.ndb" ; then echo "Genome blast database already exist" ; else makeblastdb -in $genome -dbtype nucl ; fi 
if test -f "$genome.fai" ; then echo "Genome fai file already exist" ; else samtools faidx $genome ; fi 


if [ $multiple_exon_search == "FALSE" ] ; then

	####Perform a tblastn with known V1R genes against the genome, with an evalue of 1e-5


	tblastn -query $V1R_database -db $genome -evalue $evalue -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out V1R_vs_Genome.blastn -num_threads $number_of_thread 

	cut -f1 V1R_vs_Genome.blastn  | sort | uniq -c | sed 's/^ *//g'| sed 's/ /,/g' | awk -F, ' ($1 < 3000) ' | cut -f2 -d "," > good_tblastn_id ; for i in `cat good_tblastn_id` ; do grep "^$i	" V1R_vs_Genome.blastn >> Filtered_V1R_vs_Genome.blastn ; done 
	mv Filtered_V1R_vs_Genome.blastn V1R_vs_Genome.blastn

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
	
	
	####Extract ORF that are atleast 750bp in size and rename output sequences so that the start coordinate is the start of the ORF (same for end)
	getorf -sequence Best_hits_filtered.fasta -outseq orf_list.fasta -minsize 810 -find 3
	
	
	
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
	
	
	###We apply a first filter to ORFs : blastp sequences against a database contining OR, TAAR, V2R, V1R and other GPCR genes
	
	blastp -query renamed_orf_list.fasta_uniq.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads $number_of_thread
	
	###Retain only genes that best match to TAAR
	grep -i "V1R-Receptor" blastp_result | cut -f1 | sort | uniq > V1R_from_blast.list
	xargs samtools faidx renamed_orf_list.fasta_uniq.prot < V1R_from_blast.list > Fasta_V1R_from_blast.prot
	
	
	###Second gene filter : Lets align our sequences and perform a ML tree with known TAAR genes (max 200 optimization round)
	mafft --add Fasta_V1R_from_blast.prot --keeplength $scripts_location/Database_V1R_cdhit_70_plus_T2R.aln > Putative_V1R_plus_known_V1R_plus_outgroup.prot.aln
	iqtree -s Putative_V1R_plus_known_V1R_plus_outgroup.prot.aln -st AA -nt $number_of_thread -m JTT+F+I+G4 -redo -n 200
	
	###We use a script to remove outgroup sequences that are not TAAR
	
	cp $scripts_location/V1R_sequences.id ./
	cp $scripts_location/T2R_id.txt ./
	Rscript $scripts_location/Tree_parser_v1r.R
	
	
	###Remaining sequences are TAAR. Lets remove 100% identical sequences and remove complete sequences but with ambigous nucleotides
	xargs samtools faidx renamed_orf_list.fasta_uniq.fa < Current_species_V1R.txt > Functionnal_V1Rs_multifasta_singleexon.fa
	
	#remove ambigous sequences
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_V1Rs_multifasta_singleexon.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Functionnal_V1Rs_multifasta_singleexon.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_V1Rs_multifasta_singleexon.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > unclear_Functionnal_V1Rs_multifasta_singleexon.fa


	#Translate functional genes
	transeq clear_Functionnal_V1Rs_multifasta_singleexon.fa clear_Functionnal_V1Rs_multifasta_singleexon.prot ; sed -i 's/_1$//g' clear_Functionnal_V1Rs_multifasta_singleexon.prot
	
	
	###Extract the coordinates of functionnals ORs found. Then use the script to perform the tblastn with an evalue of 1e-20
	
	grep ">" Functionnal_V1Rs_multifasta_singleexon.fa | sed 's/>//g' | sed 's/-/	/g' > Coordinates_Functionnal_V1RS.txt
	
	#if no functionnal genes found then put anything on the coordinate file
	if [ `wc -l < Coordinates_Functionnal_V1RS.txt` -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_Functionnal_V1RS.txt ; fi
	
	
	#perform second tblastn with found genes
	tblastn -query clear_Functionnal_V1Rs_multifasta_singleexon.prot -db $genome -evalue $evalue -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out tblastn_functionnal_V1R_vs_genome_current_species.tblastn -num_threads $number_of_thread

	#Remove problematic sequences
	cut -f1 tblastn_functionnal_V1R_vs_genome_current_species.tblastn  | sort | uniq -c | sed 's/^ *//g'| sed 's/ /,/g' | awk -F, ' ($1 < 4000) ' | cut -f2 -d "," > good_tblastn_id ; for i in `cat good_tblastn_id` ; do grep "$i" tblastn_functionnal_V1R_vs_genome_current_species.tblastn >> Filtered_tblastn_functionnal_V1R_vs_genome_current_species.tblastn ; done
	cat V1R_vs_Genome.blastn Filtered_tblastn_functionnal_V1R_vs_genome_current_species.tblastn > tblastn_functionnal_and_known_V1R_vs_genome.tblastn


	#cat V1R_vs_Genome.blastn tblastn_functionnal_V1R_vs_genome_current_species.tblastn > tblastn_functionnal_and_known_V1R_vs_genome.tblastn
	
	cat $V1R_database clear_Functionnal_V1Rs_multifasta_singleexon.prot > Complete_V1R_db.prot
	
	
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################## Extraction of singe exon V1R pseudogenes  ####################################################################################
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
	blastx -query Pseudogenes_regions.fa -db $blast_database -max_target_seqs 1 -outfmt '6 qseqid sseqid' -out blastx_blast_regions.tsv -num_threads $number_of_thread
	grep "V1R-Receptor" blastx_blast_regions.tsv | cut -f1 | sort | uniq > V1R_best_hits_regions.tsv
	xargs samtools faidx $genome < V1R_best_hits_regions.tsv > V1R_best_hits_regions.fa
	
	
	IFS=$'\n'
	for line in `cat V1R_best_hits_regions.tsv` ; do 
		scaffold=`echo "$line" | sed 's/:/	/g' | sed 's/-/	/g'| cut -f1`
		start=`echo "$line" | sed 's/:/	/g' | sed 's/-/	/g'| cut -f2`
		stop=`echo "$line" | sed 's/:/	/g' | sed 's/-/	/g'| cut -f3`
	
		grep "$scaffold.*$start.*$stop" Pseudo_truncated_coordinates.tsv >> Pseudo_truncated_coordinates_filtered.tsv
	done
	
	
	
	#Lets loop over found regions to characterize genes
	
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
		samtools faidx Complete_V1R_db.prot $query > Gene.prot
			 
	
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
			if grep -q "C4 Alignment" Exonerate_single_exon_pseudo ; then echo "Exonerate found something" ; else tblastn -query Gene.prot -db Genomic_region.fa -evalue 1e-02 -outfmt '6 sseq' > alternative_exo.txt ; sed -e "1i>$header_name\\" alternative_exo.txt > alternative_exo.prot ; cat alternative_exo.prot >> Frameshiftless_Pseudogenes_single_exon.prot ; fi 
		
		
		
		
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
			if grep -q "C4 Alignment" Exonerate_single_exon_pseudo ; then echo "Exonerate found something" ; else tblastn -query Gene.prot -db Genomic_region.fa -evalue 1e-02 -outfmt '6 sseq' > alternative_exo.txt ; sed -e "1i>$header_name\\" alternative_exo.txt > alternative_exo.prot ; cat alternative_exo.prot >> Frameshiftless_Pseudogenes_single_exon.prot ; fi 


		fi

	done	
	
	
	
	
	#IFS=$'\n' ; for line in `cat Pseudo_truncated_coordinates.tsv` ; do scaffold=`echo $line | cut -f1` ; start=`echo $line | cut -f2` ; stop=`echo $line | cut -f3` ; query=`echo $line | cut -f4` ; scaffold_length=`grep "^$scaffold " $genome.fai | cut -f2` ; true_start=`echo $line | cut -f5` ; true_end=`echo $line | cut -f6` ; strand=`echo $line | cut -f7` ; stop_codon_state="FALSE" ; edge_state="FALSE" ; frameshift_state="FALSE" ; samtools faidx $genome $scaffold:$start-$stop > Genomic_region.fa ; makeblastdb -in Genomic_region.fa -dbtype nucl ; samtools faidx Complete_V1R_db.prot $query > Gene.prot ; blastx -query Genomic_region.fa -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastx_result -max_target_seqs 1 -num_threads 10 ; if grep -q -i "V1R-Receptor" blastx_result ; then tblastn -query Gene.prot -db Genomic_region.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe sseq' > tblast_result_stop.tsv ; tblastn -gapopen 32767 -gapextend 32767 -query Gene.prot -db Genomic_region.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe sseq' > tblast_result_fs.tsv ; cut -f12 tblast_result_fs.tsv > Frame_detection.txt ; cut -f13 tblast_result_stop.tsv > Stop_detection.txt ; number_hsp=`cat Frame_detection.txt | wc -l`  ; number_diff_frames=`cat Frame_detection.txt | sort | uniq | wc -l` ; if [ "$number_diff_frames" -gt '1' ] ; then frameshift_state="TRUE" ; fi ; number_stops=`grep -o "\*" Stop_detection.txt | wc -l` ; if [ "$number_stops" -ge '1' ] ; then stop_codon_state="TRUE" ; fi ; if [ "$true_start" -le '100' ] ; then edge_state="TRUE" ; fi  ; diff_lengths=$((scaffold_length - true_end)) ; if [ "$diff_lengths" -le '100' ] ; then edge_state="TRUE" ; fi ; extanded_start_coord=$((start - 100)) ; extanded_end_coord=$((stop + 100)) ; consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1` ; if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi ; if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi  ; if [ $strand == "+" ] ; then samtools faidx $genome $scaffold:$true_start-$true_end > temporary_rslt.fa ; grep -v ">" temporary_rslt.fa > temporary_rslt.txt ; header_name=`echo "$scaffold-$true_start-$true_end---1_exon-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'` ; sed -e "1i>$header_name\\"  temporary_rslt.txt > temporary_rslt_renamed.fa ; sed '/^[[:space:]]*$/d' temporary_rslt_renamed.fa >> Pseudogenes_single_exon.fa ; $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 Gene.prot Genomic_region.fa > Exonerate_single_exon_pseudo ; awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' Exonerate_single_exon_pseudo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa ; transeq predicted_cds.fa predicted_cds.prot ; rm predicted_cds.fa ; sed "s/>.*/>$header_name/g" predicted_cds.prot >> Frameshiftless_Pseudogenes_single_exon.prot ; else samtools faidx $genome $scaffold:$true_start-$true_end > temporary_rslt.fa ; revseq temporary_rslt.fa temporary_rslt_rev.fa ; grep -v ">" temporary_rslt_rev.fa > temporary_rslt_rev.txt ; header_name=`echo "$scaffold-$true_start-$true_end---1_exon-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'` ; sed -e "1i>$header_name\\" temporary_rslt_rev.txt > temporary_rslt_renamed.fa ; sed '/^[[:space:]]*$/d' temporary_rslt_renamed.fa >> Pseudogenes_single_exon.fa ; $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 Gene.prot Genomic_region.fa > Exonerate_single_exon_pseudo  ; awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' Exonerate_single_exon_pseudo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa ; transeq predicted_cds.fa predicted_cds.prot ; rm predicted_cds.fa ; sed "s/>.*/>$header_name/g" predicted_cds.prot >> Frameshiftless_Pseudogenes_single_exon.prot ; fi ;fi ; done
	
	# Remove pseudogenes with ambigous nucleotides or a length less than 200nt
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_single_exon.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > Pseudogenes_single_exon_clean.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_single_exon.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > Pseudogenes_single_exon_unclean.fa 
	
	## Contrary to TAAR genes, from what I observed, V1R genes don't need a tree as further filter. 
	fasta_formatter -i Pseudogenes_single_exon_clean.fa > test.fasta ; mv test.fasta Pseudogenes_single_exon_clean.fa
	
	#Rename pseudogene fasta file
	cp Pseudogenes_single_exon_clean.fa Final_Pseudogenes.fa
	
	#Add the number of exon to functional genes
	grep ">" clear_Functionnal_V1Rs_multifasta_singleexon.fa | sed 's/>//g' > single_exon_id.txt
	sed -e 's/$/---1_exons/' single_exon_id.txt > single_exon_id_edit.txt
	paste -d "\t" single_exon_id.txt single_exon_id_edit.txt > renaming_single_exon.txt
	perl $scripts_location/rename_fasta.pl renaming_single_exon.txt clear_Functionnal_V1Rs_multifasta_singleexon.fa > temporary.fasta ; mv temporary.fasta clear_Functionnal_V1Rs_multifasta_singleexon.fa
	
	
	cp clear_Functionnal_V1Rs_multifasta_singleexon.fa Combined_Functionnal_V1R.fa
	
	
	#Check if complete genes have 7tm domain determined with phobius or TMHMM
	
	
	#First check with phobius
	transeq Combined_Functionnal_V1R.fa Combined_Functionnal_V1R.prot ; sed -i 's/_1$//g' Combined_Functionnal_V1R.prot #translate CDS
	perl $scripts_location/phobius/phobius.pl -long Combined_Functionnal_V1R.prot > Phobius_verification.txt #run phonius in long mode
	grep ">" Combined_Functionnal_V1R.prot | sed 's/>//g' > gene_id.txt #extract cds id
	
	for gene in `cat gene_id.txt` ; do nb_transm=`sed '/'"$gene"'/,/\/\//!d;/\/\//q' Phobius_verification.txt | grep "TRANSMEM" | wc -l` ; echo "$gene,$nb_transm" ; done > Gene_NbTm.tsv
	awk 'BEGIN{FS=",";OFS=","}($2>=7){print $1;}' Gene_NbTm.tsv > Phobius_genes_with_7tm.txt
	awk 'BEGIN{FS=",";OFS=","}($2<7){print $1;}' Gene_NbTm.tsv > Phobius_genes_without_7tm.txt
	
	#Now, with TMHMM
	
	$scripts_location/tmhmm-2.0c/bin/tmhmm Combined_Functionnal_V1R.prot > tmhmm_verification.txt
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
	
	perl $scripts_location/rename_fasta.pl renaming_file_tm.txt Combined_Functionnal_V1R.fa > temporary.fasta ; mv temporary.fasta Combined_Functionnal_V1R.fa
	
	
	#merge the two file containing ambigous sequences
	
	cat unclear_Functionnal_V1Rs_multifasta_singleexon.fa Pseudogenes_single_exon_unclean.fa  > Ambigous_V1R.fasta
	
	
	
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
	
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Ambigous_V1R.fasta | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Ambigous_V1R_uniq.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Final_Pseudogenes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Final_Pseudogenes_uniq.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Combined_Functionnal_V1R.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Combined_Functionnal_V1R_uniq.fa
	
	
	
	nb_seq=`grep -c ">" Ambigous_V1R_uniq.fa`
	if [ "$nb_seq" -gt "0" ] ; then
		grep ">" Ambigous_V1R_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_ambigous_final.tsv
	fi
	
	nb_seq=`grep -c ">" Final_Pseudogenes_uniq.fa`
	if [ "$nb_seq" -gt "0" ] ; then
		grep ">" Final_Pseudogenes_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_pseudogenes_final.tsv
	fi
	
	nb_seq=`grep -c ">" Combined_Functionnal_V1R_uniq.fa`
	if [ "$nb_seq" -gt "0" ] ; then
		grep ">" Combined_Functionnal_V1R_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_genes_final.tsv
	fi
	
	
	
	
	Rscript $scripts_location/Remove_redundancy.R
	
	
	
	IFS=$'\n'
	
	
	for line in `cat best_genes_functionnal.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Combined_Functionnal_V1R_uniq.fa | sed 's/>//g' >> functionnal_to_keep.txt ; done
	for line in `cat best_genes_ambigous.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Ambigous_V1R_uniq.fa | sed 's/>//g' >> ambigous_to_keep.txt ; done
	for line in `cat best_genes_pseudogenes.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Final_Pseudogenes_uniq.fa | sed 's/>//g' >> pseudogenes_to_keep.txt ; done
	
	
	xargs samtools faidx Combined_Functionnal_V1R_uniq.fa < functionnal_to_keep.txt > FINAL_Functionnal_V1R.fa
	xargs samtools faidx Ambigous_V1R_uniq.fa < ambigous_to_keep.txt > FINAL_Ambigous_V1R.fa
	xargs samtools faidx Final_Pseudogenes_uniq.fa < pseudogenes_to_keep.txt > FINAL_Pseudogenes_V1R.fa
	
	
	#Also classify genes without 7tm as pseudogenes
	
	if [ $tm_filter == "TRUE" ] ; then

		grep ">" FINAL_Functionnal_V1R.fa | grep "phobius\|tmhmm" | sed 's/>//g' > 7tm_genes
		grep ">" FINAL_Functionnal_V1R.fa | grep -v "phobius\|tmhmm" | sed 's/>//g' > non_7tm_genes
		xargs samtools faidx FINAL_Functionnal_V1R.fa < 7tm_genes > FINAL_Functionnal_V1R_7tm.fa 
		xargs samtools faidx FINAL_Functionnal_V1R.fa < non_7tm_genes >> FINAL_Pseudogenes_V1R.fa
	
	else
		cp FINAL_Functionnal_V1R.fa FINAL_Functionnal_V1R_7tm.fa
	fi





	transeq FINAL_Functionnal_V1R_7tm.fa FINAL_Functionnal_V1R_7tm.prot ; sed -i 's/_1$//g' FINAL_Functionnal_V1R_7tm.prot
	blastp -query FINAL_Functionnal_V1R_7tm.prot -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastp_result_custom -max_target_seqs 1 -num_threads $number_of_thread
	grep "V1R" blastp_result_custom | cut -f1 | sort | uniq  > good_sequences
	xargs samtools faidx FINAL_Functionnal_V1R_7tm.fa < good_sequences > temp.fa
	mv temp.fa FINAL_Functionnal_V1R_7tm.fa ; rm good_sequences ; rm *.fai
	
	
	blastx -query FINAL_Pseudogenes_V1R.fa -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastx_result_custom -max_target_seqs 1 -num_threads $number_of_thread
	grep "V1R" blastx_result_custom | cut -f1 | sort | uniq  > good_sequences_p
	xargs samtools faidx FINAL_Pseudogenes_V1R.fa < good_sequences_p > temp_p.fa
	mv temp_p.fa FINAL_Pseudogenes_V1R.fa ; rm good_sequences_p ; rm *.fai
	
	blastx -query FINAL_Ambigous_V1R.fa -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastx_result_custom -max_target_seqs 1 -num_threads $number_of_thread
	grep "V1R" blastx_result_custom | cut -f1 | sort | uniq  > good_sequences_p
	xargs samtools faidx FINAL_Ambigous_V1R.fa < good_sequences_p > temp_p.fa
	mv temp_p.fa FINAL_Ambigous_V1R.fa ; rm good_sequences_p ; rm *.fai
	
	
	#We have three final file :
	#All potentially functionnal genes : FINAL_Functionnal_V1R_7tm.fa
	#Probably pseudogenes or edge genes : FINAL_Pseudogenes_V1R.fa
	#Ambigous sequences : FINAL_Ambigous_V1R.fa
	
	nb_functionnal=`grep -c ">" FINAL_Functionnal_V1R_7tm.fa`
	nb_pseudo_edge=`grep -c ">" FINAL_Pseudogenes_V1R.fa`
	nb_ambigous=`grep -c ">" FINAL_Ambigous_V1R.fa`
	
	echo "Search of V1R is finished. There are $nb_functionnal potentially functionnal genes, $nb_pseudo_edge pseudogenes or fragments and $nb_ambigous ambigous sequences"
	
	echo "$nb_functionnal	$nb_pseudo_edge	$nb_ambigous" > Results_NbF_NbP_NbA_summary.txt
	
	
	
	dt=$(date '+%d/%m/%Y %H:%M:%S');
	echo "$dt"


fi

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## END SINGLE EXON MODE ########################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


if [ $multiple_exon_search == "TRUE" ] ; then



	#Perform tblastn using known V1R genes against the genome with an evalue of 1e-10
	
	tblastn -query $V1R_database -db $genome -evalue $evalue -outfmt 6 -out V1R_vs_Genome.blastn -num_threads $number_of_thread

	cut -f1 V1R_vs_Genome.blastn  | sort | uniq -c | sed 's/^ *//g'| sed 's/ /,/g' | awk -F, ' ($1 < 3000) ' | cut -f2 -d "," > good_tblastn_id ; for i in `cat good_tblastn_id` ; do grep "^$i	" V1R_vs_Genome.blastn >> Filtered_V1R_vs_Genome.blastn ; done 
	mv Filtered_V1R_vs_Genome.blastn V1R_vs_Genome.blastn
	
	#Lets launch a Rscript that will merge all blast hits 
	
	Rscript $scripts_location/Rscript_merge_blast_hits_7avril.R
	
	xargs samtools faidx $genome < Blast_nonoverlapping.tsv > Blast_nonoverlapping.fasta
	
	#Retain only besthit that best match to a V1R gene
	blastx -query Blast_nonoverlapping.fasta -db $blast_database -max_target_seqs 1 -outfmt '6 qseqid sseqid' -out blastx_blast_regions.tsv -num_threads $number_of_thread
	grep "V1R-Receptor" blastx_blast_regions.tsv | cut -f1 | sort | uniq > V1R_best_hits.txt
	[ -e V1R_Regions.tsv ] && rm V1R_Regions.tsv
	for i in `cat V1R_best_hits.txt` ; do grep "$i" Blast_nonoverlapping.tsv >> V1R_Regions.tsv ; done
	
	
	
	#Extend all best hits by 10000bp upstream and downstream . Result file : Potential_V1R_regions.tsv
	Rscript $scripts_location/Rscript_merge_filter_extend_blast_hit_7avril.R
	
	
	#Split the V1R database and launch exonerate with these sequences against potential V1R regions (max intron length : 30000bp)
	
	mkdir Splitted_db
	$scripts_location/exonerate-2.2.0-x86_64/bin/fastasplit -f $V1R_database -c $number_of_thread --output Splitted_db
	
	
	
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################## Search for V1R genes in a loop  #############################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	
	
	
	#re-initialize files
	
	[ -e Potential_multiple_exon_CDS.fa ] && rm Potential_multiple_exon_CDS.fa 
	[ -e Pseudogenes_multiple_exon.fa ] && rm Pseudogenes_multiple_exon.fa 
	[ -e No_V1R_genes_coordinates.txt ] && rm No_V1R_genes_coordinates.txt
	[ -e Frameshift_less_Pseudogenes.fa ] && rm Frameshift_less_Pseudogenes.fa
	
	#Start the loop to search for V1R genes
	
	current_nb_sequences=1
	previous_iteration_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`
	number_regions_blast=`grep "[0-9]" Potential_V1R_regions.tsv | wc -l`
	
	
	while [ "$current_nb_sequences" -gt "$previous_iteration_nb_sequences" ] && [ "$number_regions_blast" -gt "0" ] ; do
	
	
		previous_iteration_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`
		
		#Extract identified regions in a fasta file
	
		xargs samtools faidx $genome < Potential_V1R_regions.tsv > Potential_V1R_regions.fa
	
		mkdir Exonerate_raw_results_folder
		
		
		for i in Splitted_db/* ; do
			file_name=`echo $i | sed 's/Splitted_db\///g'`
			sbatch -W -c 4 --qos=6hours --wrap="$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Splitted_db/$file_name Potential_V1R_regions.fa > Exonerate_raw_results_folder/$file_name.exo.rslt ; sleep 10" &
		done
		
		echo "Exonerate running -- Wait"
		
		wait
		
		echo "Exonerate research done"
		
		
		#Merge exonerate results
		
		cat Exonerate_raw_results_folder/*.exo.rslt > Exonerate_results.txt
		
		
		#extract vulgar lines
		
		grep "vulgar" Exonerate_results.txt > vulgar_lines.txt 
		
		
		#extract interesting columns of vulgar lines
		#query, query_start, query_end, scaffold, scaffold_start, scaffold_end, strand, exonerate_score
		
		sed 's/vulgar: //g' vulgar_lines.txt | cut -f1,2,3,5,6,7,8,9 -d " " > vulgar_lines_parsed.txt
		
		
		#count the number of introns using vulgar lines
		
		
		IFS=$'\n'
		
		awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/ I /,"") "\t" NR}' vulgar_lines.txt > number_introns_per_line.txt
		grep -v "count" number_introns_per_line.txt | cut -f1  > intron_numbers.txt 
		
		#add the intron number to vulgar lines 
		
		paste -d " " vulgar_lines_parsed.txt intron_numbers.txt > vulgar_lines_intron_numbers.txt
		
		
		##Add informations about the best blastp results of each exonerate predicted genes
		
		
		sed -n '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/p' Exonerate_results.txt | sed 's/C4 Alignment:.*//g' | sed 's/Hostname:.*//g' | sed 's/Command line:.*//g' | sed 's/^--$//g' | sed 's/-- completed exonerate analysis.*//g' | sed 's/# --- END OF GFF DUMP ---//g' | sed 's/^#$/>seq_to_rename/g' > List_exonerate_cds.fasta #extract all predicted genes sequences
		transeq List_exonerate_cds.fasta List_exonerate_cds.prot #translate sequences
		sed 's/ /_/g' vulgar_lines_intron_numbers.txt > sequences_names.txt #extract exonerate vulgar line to rename sequences
		awk '/^>/ { printf("%s_%s\n",$0,i++);next;} { print $0;}' List_exonerate_cds.prot > List_exonerate_cds_renamed.prot #first round of rename
		grep ">" List_exonerate_cds_renamed.prot | sed 's/>//g' > old_names.txt #extract names
		paste -d "\t" old_names.txt sequences_names.txt > renaming_file #crate a file for rename_fasta.pl
		perl $scripts_location/rename_fasta.pl renaming_file List_exonerate_cds_renamed.prot > List_exonerate_cds.prot #completely rename sequences with the exonerate vulgar line
		
		#Perform the blastp
		blastp -query List_exonerate_cds.prot -db $scripts_location/Database_V1R.prot -outfmt '6 qseqid sseqid evalue' -out all_blastp.txt -max_target_seqs 1 -num_threads $number_of_thread
		
		#Extract the information
		[ -e all_blastp_parsed.txt ] && rm all_blastp_parsed.txt
		for i in `cat sequences_names.txt` ; do if grep -q "$i" all_blastp.txt ; then grep -m1 "$i" all_blastp.txt | cut -f2,3 >> all_blastp_parsed.txt ; else echo "NoQuery	99999" >> all_blastp_parsed.txt ; fi ; done 
		sed -i 's/	/ /g' all_blastp_parsed.txt
		paste -d " " vulgar_lines_intron_numbers.txt all_blastp_parsed.txt > vulgar_lines_intron_numbers_blastrslt.txt
		
		
		
		#Parse exonerate results. Find the best exonerate results that are most likely complete genes or pseudoogenes, and not overlapping
		
		Rscript $scripts_location/Parse_exonerate_results.R #result file : Parsed_exonerate_gene_regions.tsv
		
		
	
	
	
		nb_row_parsed_exonerate=`wc -l < Parsed_exonerate_gene_regions.tsv`
		if [ "$nb_row_parsed_exonerate" -gt "0" ] ; then
	
	
			IFS=$'\n'
			
			for line in `cat Parsed_exonerate_gene_regions.tsv` ; do
				
				query=`echo "$line" | cut -f7`
				scaffold=`echo "$line" | cut -f1`
				scaff_start=`echo "$line" | cut -f2`
				scaff_end=`echo "$line" | cut -f3`
			
				echo "$scaffold:$scaff_start-$scaff_end	$query"
			
			done > Correct_coordinates_for_exonerate.tsv
			
			
			#Let's now predict genes on these regions !
			
			IFS=$'\n'
			
			
			mkdir Genes_predictions
			
			
			for line in `cat Correct_coordinates_for_exonerate.tsv` ; do 
				
				scaffold_s_e=`echo "$line" | cut -f1`
				best_query=`echo "$line" | cut -f2`
				scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'`
			
			
				samtools faidx $genome $scaffold_s_e > scaffold.fa
				sed -i 's/:/-/g' scaffold.fa
			
				samtools faidx $scripts_location/Database_V1R.prot $best_query > query.prot
			
				$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
				cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta
			
			
				#extract only the best result if there are two with the same score
			
				if [ `grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate` -ge 2 ] ; then
					sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos
					sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence
					cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate
				fi
			
			
			
			done
			
			
			### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
			
			
			#Result folder
			mkdir Filtered_predictions	
			
			for file in Genes_predictions/*.exonerate ; do
			
				#extract some infos from file name
				file_name=`echo "$file" | sed 's/.*\///g'`
				file_name_reduced=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
				fasta_file_name=`echo "$file_name" | sed 's/exonerate/fasta/g'`
				initial_header=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g'`
			
				#Test if the predicted gene is a V1R gene or not
				awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
				transeq predicted_cds.fa predicted_cds.prot
				blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 10
				
			
				#Lets continue only if the best match is an OR
				#if grep -q -i "olfactory\|odorant" blastp_result ; then 
				if grep -q -i "V1R-Receptor" blastp_result ; then
			
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
					
					#file=Genes_predictions/NC_019879.2-28438421-28440954.exonerate 
			
			
						#If strand is minus, then the first position is:
						target_end=$((first_hit_range + 1))
						#And we will went to extend this by 500bp to be sure to have the potentiel start codon
						target_extanded_end=$((first_hit_range + 500))
					
						#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
						samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
						samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
					
						#remove fasta header of extanded region files
						grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
						grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
					
						#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
			
						grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
						for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
						grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
			
					
						#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
						cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
						sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
						getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 810 -find 3
						if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then sequence_to_grep=`grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/>//g' | sed 's/ .*//g'`  ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $sequence_to_grep > temporary ; mv temporary Filtered_predictions/$file_name_reduced.ORF ; rm Filtered_predictions/$file_name_reduced.ORF.fai ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
					
			
			
						#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
						if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
					
							transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
							$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
							if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi
	
	
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
							number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | wc -l`
							if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
			
			
			
			
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
							query_total_length=`grep -m1 "$query_name" $scripts_location/Database_V1R.prot.fai | cut -f2`
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
					
							extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
					
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
					
								samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
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
							number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | wc -l`
							if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
							if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
		
					
						fi
			
			
			
					#Lets make the same steps with slight modifications for the + strand
					elif [ $strand == "+" ] ; then 
					
						#If strand is minus, then the first position is:
						target_end=$((second_hit_range + 1))
						#And we will went to extend this by 500bp to be sure to have the potentiel start codon
						target_extanded_end=$((second_hit_range + 500))
					
						#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
						samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
						samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
					
						#remove fasta header of extanded region files
						grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
						grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
					
						#Extract the CDS sequence predicted by exonerate and remove fasta header
			
						grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
						for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
						grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
					
						#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
						cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
						sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
						getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 810 -find 3 -reverse FALSE
					
					
						#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
						if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
					
							transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
							$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
							if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi
	
	
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
							number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | wc -l`
							if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
			
			
			
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
							query_total_length=`grep -m1 "$query_name" $scripts_location/Database_V1R.prot.fai | cut -f2`
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
					
							extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
					
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
					
								samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
					
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
							number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | wc -l`
							if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
							if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
		
		
						fi
					fi
			
				else echo "$initial_header" >> No_V1R_genes_coordinates.txt
			
				fi
			
			done
		fi
		
		
		#Now that we have filtered all our results, we can concatenate the results
		
		for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
		cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
		cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
		cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa
	 
		#Extract coordinates of found genes
		grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_already_examined.tsv
		grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
		sed 's/-/	/g' No_V1R_genes_coordinates.txt >> Coordinates_already_examined.tsv
		if [ `wc -l < Coordinates_already_examined.tsv` -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_already_examined.tsv ; fi
		
		current_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`
		
		#re-process blast result to find potential V1R regions exclusing already found genes
		
		Rscript $scripts_location/Rscript_merge_filter_extend_blast_hit_Second_7avril.R
	
	
	
		rm -r Filtered_predictions/
		rm -r Genes_predictions/
		rm -r Exonerate_raw_results_folder/
		rm Parsed_exonerate_gene_regions.tsv
		
	
		if test -f "Potential_V1R_regions.tsv" ; then number_regions_blast=`grep "[0-9]" Potential_V1R_regions.tsv | wc -l` ; else number_regions_blast=0 ; fi
	
		
	done
	
	
	rm -r Filtered_predictions/
	rm -r Genes_predictions/
	rm Parsed_exonerate_gene_regions.tsv
	
	
	
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################## Search for remaining V1R genes with size a bit below  #######################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	
	
	Rscript $scripts_location/Parse_exonerate_results_second_3Avril.R
	
	
	nb_row_parsed_exonerate=`wc -l < Parsed_exonerate_gene_regions.tsv`
	
	
	if [ "$nb_row_parsed_exonerate" -gt "0" ] ; then
	
		IFS=$'\n'
		
		for line in `cat Parsed_exonerate_gene_regions.tsv` ; do
			
			query=`echo "$line" | cut -f7`
			scaffold=`echo "$line" | cut -f1`
			scaff_start=`echo "$line" | cut -f2`
			scaff_end=`echo "$line" | cut -f3`
		
			echo "$scaffold:$scaff_start-$scaff_end	$query"
		
		done > Correct_coordinates_for_exonerate.tsv
		
		
		#Let's now predict genes on these regions !
		
		IFS=$'\n'
		
		
		mkdir Genes_predictions
		
		
		for line in `cat Correct_coordinates_for_exonerate.tsv` ; do 
			
			scaffold_s_e=`echo "$line" | cut -f1`
			best_query=`echo "$line" | cut -f2`
			scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'`
		
		
			samtools faidx $genome $scaffold_s_e > scaffold.fa
			sed -i 's/:/-/g' scaffold.fa
		
			samtools faidx $scripts_location/Database_V1R.prot $best_query > query.prot
		
			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
			cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta
		
		
			#extract only the best result if there are two with the same score
		
			if [ `grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate` -ge 2 ] ; then
				sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos
				sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence
				cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate
			fi
		
		
		
		done
		
		
		### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
		
		
		#Result folder
		mkdir Filtered_predictions	
		
		for file in Genes_predictions/*.exonerate ; do
		
			#extract some infos from file name
			file_name=`echo "$file" | sed 's/.*\///g'`
			file_name_reduced=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
			fasta_file_name=`echo "$file_name" | sed 's/exonerate/fasta/g'`
			initial_header=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g'`
		
			#Test if the predicted gene is a V1R gene or not
			awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
			transeq predicted_cds.fa predicted_cds.prot
			blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 10
			
		
			#Lets continue only if the best match is an OR
			#if grep -q -i "olfactory\|odorant" blastp_result ; then 
			if grep -q -i "V1R-Receptor" blastp_result ; then
		
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
				
				#file=Genes_predictions/NC_019879.2-28438421-28440954.exonerate 
		
		
					#If strand is minus, then the first position is:
					target_end=$((first_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((first_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
				
					#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
		
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
		
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 810 -find 3
					if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then sequence_to_grep=`grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/>//g' | sed 's/ .*//g'`  ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $sequence_to_grep > temporary ; mv temporary Filtered_predictions/$file_name_reduced.ORF ; rm Filtered_predictions/$file_name_reduced.ORF.fai ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
				
		
		
					#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
					if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
						$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi
	
	
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
						number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | wc -l`
						if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
		
		
		
		
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
						query_total_length=`grep -m1 "$query_name" $scripts_location/Database_V1R.prot.fai | cut -f2`
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
				
						extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
				
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
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
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
						number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | wc -l`
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi
		
		
				
					fi
		
		
		
				#Lets make the same steps with slight modifications for the + strand
				elif [ $strand == "+" ] ; then 
				
					#If strand is minus, then the first position is:
					target_end=$((second_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((second_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
				
					#Extract the CDS sequence predicted by exonerate and remove fasta header
		
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 810 -find 3 -reverse FALSE
				
				
					#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
					if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
						$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi
	
	
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
						number_blast_hit=`tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | wc -l`
						if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
		
		
		
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
						query_total_length=`grep -m1 "$query_name" $scripts_location/Database_V1R.prot.fai | cut -f2`
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
				
						extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
				
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
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
				
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
						number_blast_hit=`tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue 1e-10 -outfmt 6 | wc -l`
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi
		
		
					fi
				fi
		
			else echo "$initial_header" >> No_V1R_genes_coordinates.txt
		
			fi
		
		done
	
	fi
	
	#Now that we have filtered all our results, we can concatenate the results
	
	for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
	cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
	cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
	cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa
	
	#Extract coordinates of found genes
	grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_already_examined.tsv
	grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
	sed 's/-/	/g' No_V1R_genes_coordinates.txt >> Coordinates_already_examined.tsv
	
	current_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`
	
	#re-process blast result to find potential V1R regions exclusing already found genes
	
	rm -r Filtered_predictions/
	rm -r Genes_predictions/
	rm Parsed_exonerate_gene_regions.tsv
	
	
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################## Search for V1R pseudogenes  #################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	
	Rscript $scripts_location/Parse_exonerate_results_third_3Avril.R
	
	
	
	nb_row_parsed_exonerate=`wc -l < Parsed_exonerate_gene_regions.tsv`
	
	
	if [ "$nb_row_parsed_exonerate" -gt "0" ] ; then
	
		IFS=$'\n'
		
		for line in `cat Parsed_exonerate_gene_regions.tsv` ; do
			
			query=`echo "$line" | cut -f7`
			scaffold=`echo "$line" | cut -f1`
			scaff_start=`echo "$line" | cut -f2`
			scaff_end=`echo "$line" | cut -f3`
		
			echo "$scaffold:$scaff_start-$scaff_end	$query"
		
		done > Correct_coordinates_for_exonerate.tsv
		
		
		#Let's now predict genes on these regions !
		
		IFS=$'\n'
		
		
		mkdir Genes_predictions
		
		
		for line in `cat Correct_coordinates_for_exonerate.tsv` ; do 
			
			scaffold_s_e=`echo "$line" | cut -f1`
			best_query=`echo "$line" | cut -f2`
			scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'`
		
		
			samtools faidx $genome $scaffold_s_e > scaffold.fa
			sed -i 's/:/-/g' scaffold.fa
		
			samtools faidx $scripts_location/Database_V1R.prot $best_query > query.prot
		
			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
			cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta
		
		
			#extract only the best result if there are two with the same score
		
			if [ `grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate` -ge 2 ] ; then
				sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos
				sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence
				cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate
			fi
		
		
		
		done
		
		
		### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
		
		
		#Result folder
		mkdir Filtered_predictions	
		
		for file in Genes_predictions/*.exonerate ; do
		
			#extract some infos from file name
			file_name=`echo "$file" | sed 's/.*\///g'`
			file_name_reduced=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
			fasta_file_name=`echo "$file_name" | sed 's/exonerate/fasta/g'`
			initial_header=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g'`
		
			#Test if the predicted gene is a V1R gene or not
			awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
			transeq predicted_cds.fa predicted_cds.prot
			blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 10
			
		
			#Lets continue only if the best match is an OR
			#if grep -q -i "olfactory\|odorant" blastp_result ; then 
			if grep -q -i "V1R-Receptor" blastp_result ; then
		
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
				
				#file=Genes_predictions/NC_019879.2-28438421-28440954.exonerate 
		
		
					#If strand is minus, then the first position is:
					target_end=$((first_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((first_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
				
					#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
		
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
		
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 810 -find 3
					if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then sequence_to_grep=`grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/>//g' | sed 's/ .*//g'`  ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $sequence_to_grep > temporary ; mv temporary Filtered_predictions/$file_name_reduced.ORF ; rm Filtered_predictions/$file_name_reduced.ORF.fai ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
				
		
		
					#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
					if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
						$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi
	
	
						extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
						cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
						cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=`grep "	exon	" verif_coord.exo | wc -l`
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
				
		
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
						query_total_length=`grep -m1 "$query_name" $scripts_location/Database_V1R.prot.fai | cut -f2`
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
				
						extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
				
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
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
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
		
		
				
					fi
		
		
		
				#Lets make the same steps with slight modifications for the + strand
				elif [ $strand == "+" ] ; then 
				
					#If strand is minus, then the first position is:
					target_end=$((second_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((second_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
				
					#Extract the CDS sequence predicted by exonerate and remove fasta header
		
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 810 -find 3 -reverse FALSE
				
				
					#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
					if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
						$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi
	
	
						extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
						cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
						cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=`grep "	exon	" verif_coord.exo | wc -l`
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
		
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
						query_total_length=`grep -m1 "$query_name" $scripts_location/Database_V1R.prot.fai | cut -f2`
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
				
						extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
				
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
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
				
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
		
					fi
				fi
		
			else echo "$initial_header" >> No_V1R_genes_coordinates.txt
		
			fi
		
		done
	fi
	
	
	#Now that we have filtered all our results, we can concatenate the results
	
	for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
	cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
	cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
	cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa
	
	#Extract coordinates of found genes
	grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_already_examined.tsv
	grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
	sed 's/-/	/g' No_V1R_genes_coordinates.txt >> Coordinates_already_examined.tsv
	
	current_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`
	
	#re-process blast result to find potential V1R regions exclusing already found genes
	
	rm -r Filtered_predictions/
	rm -r Genes_predictions/
	rm Parsed_exonerate_gene_regions.tsv
	
	
	
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################## Final round to find V1R pseudogenes with length below previous iteration  ###################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	
	
	Rscript $scripts_location/Parse_exonerate_results_final_3Avril.R
	
	nb_row_parsed_exonerate=`wc -l < Parsed_exonerate_gene_regions.tsv`
	
	
	if [ "$nb_row_parsed_exonerate" -gt "0" ] ; then	
	
		IFS=$'\n'
		
		for line in `cat Parsed_exonerate_gene_regions.tsv` ; do
			
			query=`echo "$line" | cut -f7`
			scaffold=`echo "$line" | cut -f1`
			scaff_start=`echo "$line" | cut -f2`
			scaff_end=`echo "$line" | cut -f3`
		
			echo "$scaffold:$scaff_start-$scaff_end	$query"
		
		done > Correct_coordinates_for_exonerate.tsv
		
		
		#Let's now predict genes on these regions !
		
		IFS=$'\n'
		
		
		mkdir Genes_predictions
		
		
		for line in `cat Correct_coordinates_for_exonerate.tsv` ; do 
			
			scaffold_s_e=`echo "$line" | cut -f1`
			best_query=`echo "$line" | cut -f2`
			scaffold_s_e_n=`echo "$line" | cut -f1 | sed 's/:/-/g'`
		
		
			samtools faidx $genome $scaffold_s_e > scaffold.fa
			sed -i 's/:/-/g' scaffold.fa
		
			samtools faidx $scripts_location/Database_V1R.prot $best_query > query.prot
		
			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
			cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta
		
		
			#extract only the best result if there are two with the same score
		
			if [ `grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate` -ge 2 ] ; then
				sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos
				sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence
				cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate
			fi
		
		
		
		done
		
		
		
		### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
		
		
		#Result folder
		mkdir Filtered_predictions	
		
		for file in Genes_predictions/*.exonerate ; do
		
			#extract some infos from file name
			file_name=`echo "$file" | sed 's/.*\///g'`
			file_name_reduced=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g'`
			fasta_file_name=`echo "$file_name" | sed 's/exonerate/fasta/g'`
			initial_header=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g'`
		
			#Test if the predicted gene is a V1R gene or not
			awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
			transeq predicted_cds.fa predicted_cds.prot
			blastp -query predicted_cds.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 10
			
		
			#Lets continue only if the best match is an OR
			#if grep -q -i "olfactory\|odorant" blastp_result ; then 
			if grep -q -i "V1R-Receptor" blastp_result ; then
		
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
				
				#file=Genes_predictions/NC_019879.2-28438421-28440954.exonerate 
		
		
					#If strand is minus, then the first position is:
					target_end=$((first_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((first_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
				
					#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
		
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
		
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 810 -find 3
					if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then sequence_to_grep=`grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/>//g' | sed 's/ .*//g'`  ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $sequence_to_grep > temporary ; mv temporary Filtered_predictions/$file_name_reduced.ORF ; rm Filtered_predictions/$file_name_reduced.ORF.fai ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
				
		
		
					#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
					if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
						$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi
	
	
						extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
						cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
						cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=`grep "	exon	" verif_coord.exo | wc -l`
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
				
		
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
						query_total_length=`grep -m1 "$query_name" $scripts_location/Database_V1R.prot.fai | cut -f2`
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
				
						extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
				
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
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
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
		
				
					fi
		
		
		
				#Lets make the same steps with slight modifications for the + strand
				elif [ $strand == "+" ] ; then 
				
					#If strand is minus, then the first position is:
					target_end=$((second_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((second_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt
				
					#Extract the CDS sequence predicted by exonerate and remove fasta header
		
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize 810 -find 3 -reverse FALSE
				
				
					#Rename the fasta file (might also be usefull to generate a gff3 file using exonerate ? )
					if [ `grep -c ">" Filtered_predictions/$file_name_reduced.ORF` -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP
						$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi
	
	
						extracted_scaffold_start=`echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2`
						cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
						cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=`grep "	exon	" verif_coord.exo | wc -l`
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
					
		
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
						query_total_length=`grep -m1 "$query_name" $scripts_location/Database_V1R.prot.fai | cut -f2`
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
				
						extracted_scaffold_start=`grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2`
				
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
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
				
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
		
					fi
				fi
		
			else echo "$initial_header" >> No_V1R_genes_coordinates.txt
		
			fi
		
		done
	fi
	
	
	
	#Now that we have filtered all our results, we can concatenate the results
	
	for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
	cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
	cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
	cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa
	
	#Extract coordinates of found genes
	grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_already_examined.tsv
	grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
	sed 's/-/	/g' No_V1R_genes_coordinates.txt >> Coordinates_already_examined.tsv
	
	current_nb_sequences=`if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi`
	
	#re-process blast result to find potential V1R regions exclusing already found genes
	
	rm -r Filtered_predictions/
	rm -r Genes_predictions/
	rm Parsed_exonerate_gene_regions.tsv
	
	
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	##################For V1R genes , we can complete exonerate predictions with classic ORF finding  ##############################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	
	
	####Extract ORF that are atleast 810bp in size and rename output sequences so that the start coordinate is the start of the ORF (same for end)
	sed -i 's/:/-/g' Potential_V1R_regions.fa
	getorf -sequence Potential_V1R_regions.fa -outseq orf_list.fasta -minsize 810 -find 3
	
	
	nb_orf_found=`grep -c ">" orf_list.fasta`
	
	
	if [ "$nb_orf_found" -gt "0" ] ; then
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
		
		sed -i 's/$/---1_exons/g' newfastaheaders
		
		paste -d "\t" oldfastaheaders newfastaheaders > renaming_file
		
		perl $scripts_location/rename_fasta.pl renaming_file orf_list.fasta > renamed_orf_list.fasta
		
		
		###Lets translate DNA to Prot
		
		transeq renamed_orf_list.fasta renamed_orf_list.fasta_uniq.prot ; sed -i 's/_1$//g' renamed_orf_list.fasta_uniq.prot
		
		
		###We apply a blastp filter to ORFs 
		
		blastp -query renamed_orf_list.fasta_uniq.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads $number_of_thread
		
		if grep -q "V1R-Receptor" blastp_result ; then
	
			###Retain only genes that best match to V1R
			grep -i "V1R-Receptor" blastp_result | cut -f1 | sort | uniq > V1R_from_blast.list
			xargs samtools faidx renamed_orf_list.fasta < V1R_from_blast.list > renamed_orf_list_blastp.fasta
			
			
			#Concatenate with previous results
			
			cat Potential_multiple_exon_CDS.fa renamed_orf_list_blastp.fasta > Potential_V1R_completed.fa
			
			
			###Due to their potential finding in previous steps, remove identical sequences
			
			awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Potential_V1R_completed.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Potential_multiple_exon_CDS.fa
		fi
	
	fi
	
	
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################## Filter genes and pseudogenes with a blast and a phylogenetic tree  ######################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	
	#First, remove identical headers for the following steps
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_multiple_exon.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Pseudogenes_multiple_exon_uniq.fa ; mv Pseudogenes_multiple_exon_uniq.fa Pseudogenes_multiple_exon.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Frameshift_less_Pseudogenes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Frameshift_less_Pseudogenes_uniq.fa ; mv Frameshift_less_Pseudogenes_uniq.fa Frameshift_less_Pseudogenes.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Potential_multiple_exon_CDS.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Potential_multiple_exon_CDS_uniq.fa ; mv Potential_multiple_exon_CDS_uniq.fa Potential_multiple_exon_CDS.fa
	
	
	# We now have all predicted genes in the files Pseudogenes_multiple_exon.fa and Potential_multiple_exon_CDS.fa
	

	#First , lets verify our genes via a blast against GPCR database

	#Verify that these gene fragments correspong to T2R genes with a blastx
	#transeq Potential_multiple_exon_CDS.fa Potential_multiple_exon_CDS.prot ; sed -i 's/_1$//g' Potential_multiple_exon_CDS.prot
	#blastp -query Potential_multiple_exon_CDS.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid" -out verif_blastp_result -max_target_seqs 1 -num_threads $number_of_thread
	#grep -i "V1R-Receptor" verif_blastp_result | cut -f1 | sort | uniq > V1R_from_blast.list
	#xargs samtools faidx Potential_multiple_exon_CDS.fa < V1R_from_blast.list > Fasta_V1R_from_blast.fa ; mv Fasta_V1R_from_blast.fa Potential_multiple_exon_CDS.fa
#
	#fasta_formatter -i Pseudogenes_multiple_exon.fa  > Pseudogenes_multiple_exon_reformat.fa ; mv Pseudogenes_multiple_exon_reformat.fa Pseudogenes_multiple_exon.fa ; rm Pseudogenes_multiple_exon.fa.fai
	#blastx -query Pseudogenes_multiple_exon.fa -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid" -out blastx_result -max_target_seqs 1 -num_threads $number_of_thread
	#grep -i "V1R-Receptorr" blastx_result | cut -f1 | sort | uniq > Pseudo_V1R_from_blast.list
	#xargs samtools faidx Pseudogenes_multiple_exon.fa < Pseudo_V1R_from_blast.list > Pseudogenes_multiple_exon_blastx.fa ; mv Pseudogenes_multiple_exon_blastx.fa Pseudogenes_multiple_exon_reformat.fa

	

	#Then lets verify our genes with a phylogenetic tree rooted with T1R and CasR genes
	
	
	grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' > Unverified_pseudo_id.txt 
	fasta_formatter -i Pseudogenes_multiple_exon.fa  > Pseudogenes_multiple_exon_reformat.fa
	
	#remove sequences with less than 120aa, it does not permit to discriminate below...
	#seqkit seq -m 120 Frameshift_less_Pseudogenes.fa > Frameshift_less_Pseudogenes_length.fa
	cp Frameshift_less_Pseudogenes.fa Frameshift_less_Pseudogenes_length.fa
	
	nb_seq=`grep -c ">" Frameshift_less_Pseudogenes_length.fa`
	if [ "$nb_seq" -gt "0" ] ; then
		mafft --add Frameshift_less_Pseudogenes_length.fa --keeplength $scripts_location/Database_V1R_cdhit_70_plus_T2R.aln > ALL_verification_alignment.aln
	else
		cp $scripts_location/Database_V1R_cdhit_70_plus_T2R.aln ./ALL_verification_alignment.aln
	fi
	
	
	#also add functionnal genes
	
	
	
	nb_seq=`grep -c ">" Potential_multiple_exon_CDS.fa`
	if [ "$nb_seq" -gt "0" ] ; then
		transeq Potential_multiple_exon_CDS.fa Potential_multiple_exon_CDS.prot ; sed -i 's/_1$//g' Potential_multiple_exon_CDS.prot
		mafft --add Potential_multiple_exon_CDS.prot --keeplength ALL_verification_alignment.aln > Final_ALL_verification_alignment.aln
	else
		cp ALL_verification_alignment.aln Final_ALL_verification_alignment.aln
	fi
	
	
	#lets start a tree
	
	iqtree -s Final_ALL_verification_alignment.aln -st AA -nt $number_of_thread -m JTT+F+I+G4 -redo -n 200
	
	#parse the tree
	
	cp $scripts_location/V1R_sequences.id ./
	cp $scripts_location/T2R_id.txt ./
	Rscript $scripts_location/Tree_parser.R #This script keep only genes clustering with known V1R genes : Current_species_V1R.txt
	
	
	grep "FALSE\|TRUE" Current_species_V1R.txt  > V1R_pseudos.id
	grep -v "FALSE\|TRUE" Current_species_V1R.txt  > V1R_functionnal.id
	
	
	xargs samtools faidx Potential_multiple_exon_CDS.fa < V1R_functionnal.id > Functionnal_V1R_genes.fa
	xargs samtools faidx Pseudogenes_multiple_exon_reformat.fa < V1R_pseudos.id > Pseudogenes_V1R_genes.fa
	
	
	## Now lets filter these files! 
	#clear ambigous sequences
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_V1R_genes.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Functionnal_V1R_genes.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_V1R_genes.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Pseudogenes_V1R_genes.fa
	
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_V1R_genes.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > unclear_Functionnal_V1R_genes.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_V1R_genes.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > unclear_Pseudogenes_V1R_genes.fa
	
	cat unclear_Functionnal_V1R_genes.fa unclear_Pseudogenes_V1R_genes.fa > FINAL_Ambigous_V1R.fa


	#Check 7tm domains
	
	
	#Check if complete genes have 7tm domain determined with phobius or TMHMM
	
	
	#First check with phobius
	transeq clear_Functionnal_V1R_genes.fa clear_Functionnal_V1R_genes.prot ; sed -i 's/_1$//g' clear_Functionnal_V1R_genes.prot #translate CDS
	perl $scripts_location/phobius/phobius.pl -long clear_Functionnal_V1R_genes.prot > Phobius_verification.txt #run phonius in long mode
	grep ">" clear_Functionnal_V1R_genes.prot | sed 's/>//g' > gene_id.txt #extract cds id
	
	for gene in `cat gene_id.txt` ; do nb_transm=`sed '/'"$gene"'/,/\/\//!d;/\/\//q' Phobius_verification.txt | grep "TRANSMEM" | wc -l` ; echo "$gene,$nb_transm" ; done > Gene_NbTm.tsv
	awk 'BEGIN{FS=",";OFS=","}($2>=7){print $1;}' Gene_NbTm.tsv > Phobius_genes_with_7tm.txt
	awk 'BEGIN{FS=",";OFS=","}($2<7){print $1;}' Gene_NbTm.tsv > Phobius_genes_without_7tm.txt
	
	#Now, with TMHMM
	
	$scripts_location/tmhmm-2.0c/bin/tmhmm clear_Functionnal_V1R_genes.prot > tmhmm_verification.txt
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
	
	perl $scripts_location/rename_fasta.pl renaming_file_tm.txt clear_Functionnal_V1R_genes.fa > temporary.fasta ; mv temporary.fasta clear_Functionnal_V1R_genes.fa
	
	
	# If two genes/pseudogenes are overlapping due to the multiple extensions, then keep only the longest found gene/pseudogene
	
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' FINAL_Ambigous_V1R.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > FINAL_Ambigous_V1R_uniq.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' clear_Pseudogenes_V1R_genes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > clear_Pseudogenes_V1R_genes_uniq.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' clear_Functionnal_V1R_genes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > clear_Functionnal_V1R_genes_uniq.fa
	
	
	
	nb_seq=`grep -c ">" FINAL_Ambigous_V1R_uniq.fa`
	if [ "$nb_seq" -gt "0" ] ; then
		grep ">" FINAL_Ambigous_V1R_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_ambigous_final.tsv
	fi
	
	nb_seq=`grep -c ">" clear_Pseudogenes_V1R_genes_uniq.fa`
	if [ "$nb_seq" -gt "0" ] ; then
		grep ">" clear_Pseudogenes_V1R_genes_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_pseudogenes_final.tsv
	fi
	
	nb_seq=`grep -c ">" clear_Functionnal_V1R_genes_uniq.fa`
	if [ "$nb_seq" -gt "0" ] ; then
		grep ">" clear_Functionnal_V1R_genes_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_genes_final.tsv
	fi
	
	
	
	
	Rscript $scripts_location/Remove_redundancy.R


	## For teleost V1R genes, we dont check 7tm, as there is one V1R among the six canonical, on which 
	## phobius and tmhmm almost always dont predict a correct 7tm ..;
	
	
	
	IFS=$'\n'
	
	
	for line in `cat best_genes_functionnal.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" clear_Functionnal_V1R_genes_uniq.fa | sed 's/>//g' >> functionnal_to_keep.txt ; done
	for line in `cat best_genes_ambigous.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" FINAL_Ambigous_V1R_uniq.fa | sed 's/>//g' >> ambigous_to_keep.txt ; done
	for line in `cat best_genes_pseudogenes.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" clear_Pseudogenes_V1R_genes_uniq.fa | sed 's/>//g' >> pseudogenes_to_keep.txt ; done
	
	
	xargs samtools faidx clear_Functionnal_V1R_genes_uniq.fa < functionnal_to_keep.txt > FINAL_Functionnal_V1R_7tm.fa
	xargs samtools faidx FINAL_Ambigous_V1R_uniq.fa < ambigous_to_keep.txt > FINAL_Ambigous_V1R.fa
	xargs samtools faidx clear_Pseudogenes_V1R_genes_uniq.fa < pseudogenes_to_keep.txt > FINAL_Pseudogenes_V1R.fa


	transeq FINAL_Functionnal_V1R_7tm.fa FINAL_Functionnal_V1R_7tm.prot ; sed -i 's/_1$//g' FINAL_Functionnal_V1R_7tm.prot
	blastp -query FINAL_Functionnal_V1R_7tm.prot -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastp_result_custom -max_target_seqs 1 -num_threads $number_of_thread
	grep "V1R" blastp_result_custom | cut -f1 | sort | uniq  > good_sequences
	xargs samtools faidx FINAL_Functionnal_V1R_7tm.fa < good_sequences > temp.fa
	mv temp.fa FINAL_Functionnal_V1R_7tm.fa ; rm good_sequences ; rm *.fai
	
	
	blastx -query FINAL_Pseudogenes_V1R.fa -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastx_result_custom -max_target_seqs 1 -num_threads $number_of_thread
	grep "V1R" blastx_result_custom | cut -f1 | sort | uniq  > good_sequences_p
	xargs samtools faidx FINAL_Pseudogenes_V1R.fa < good_sequences_p > temp_p.fa
	mv temp_p.fa FINAL_Pseudogenes_V1R.fa ; rm good_sequences_p ; rm *.fai
	
	blastx -query FINAL_Ambigous_V1R.fa -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastx_result_custom -max_target_seqs 1 -num_threads $number_of_thread
	grep "V1R" blastx_result_custom | cut -f1 | sort | uniq  > good_sequences_p
	xargs samtools faidx FINAL_Ambigous_V1R.fa < good_sequences_p > temp_p.fa
	mv temp_p.fa FINAL_Ambigous_V1R.fa ; rm good_sequences_p ; rm *.fai
	
	
	dt=$(date '+%d/%m/%Y %H:%M:%S');
	echo "$dt"
	
	
	
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################## END MULTIPLE EXON MODE  #####################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
	################################################################################################################################################
fi

