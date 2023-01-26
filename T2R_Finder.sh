#!/bin/bash


#SBATCH --job-name=T2R   # Job name


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
T2R_database=$2
blast_database=$3
scripts_folder_location=$4 ; scripts_location=`echo "$scripts_folder_location" | sed 's/\/$//'`
maximum_intron_length=$5
number_of_thread=$6
evalue=$7
tm_filter=$8



###Makeblastdb so we can blast genes against the genome


if test -f "$genome.ndb" ; then echo "Genome blast database already exist" ; else makeblastdb -in $genome -dbtype nucl ; fi 
if test -f "$genome.fai" ; then echo "Genome fai file already exist" ; else samtools faidx $genome ; fi 

#makeblastdb -in $genome -dbtype nucl ; samtools faidx $genome


####Perform a tblastn with known T2R genes against the genome, with an evalue of 1e-5


tblastn -query $T2R_database -db $genome -evalue $evalue -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out T2R_vs_Genome.blastn -num_threads $number_of_thread 
cut -f1 T2R_vs_Genome.blastn  | sort | uniq -c | sed 's/^ *//g'| sed 's/ /,/g' | awk -F, ' ($1 < 3000) ' | cut -f2 -d "," > good_tblastn_id ; for i in `cat good_tblastn_id` ; do grep "^$i	" T2R_vs_Genome.blastn >> Filtered_T2R_vs_Genome.blastn ; done 
mv Filtered_T2R_vs_Genome.blastn T2R_vs_Genome.blastn


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Extraction of single exon T2R genes  #########################################################################################
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
getorf -sequence Best_hits_filtered.fasta -outseq orf_list.fasta -minsize 750 -find 3



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


###We apply a first filter to ORFs : blastp sequences against a database contining OR, T2R, V2R, V1R and other GPCR genes

blastp -query renamed_orf_list.fasta_uniq.prot -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads $number_of_thread

###Retain only genes that best match to T2R
grep -i "T2R-Receptor\|Taste-receptor" blastp_result | cut -f1 | sort | uniq > T2R_from_blast.list
xargs samtools faidx renamed_orf_list.fasta_uniq.prot < T2R_from_blast.list > Fasta_T2R_from_blast.prot


###Second gene filter : Lets align our sequences and perform a ML tree with known T2R genes (max 200 optimization round)
mafft --add Fasta_T2R_from_blast.prot --keeplength $scripts_location/T2R_plus_outgroup.aln > Putative_t2r_plus_known_t2r_plus_outgroup.prot.aln
iqtree -s Putative_t2r_plus_known_t2r_plus_outgroup.prot.aln -st AA -nt $number_of_thread -m JTT+F+I+G4 -redo -n 200

###We use a script to remove outgroup sequences that are not T2R
cp $scripts_location/T2R_genes_tree.id ./
cp $scripts_location/V1R_genes_tree.id ./
Rscript $scripts_location/Tree_parser.R 

xargs samtools faidx renamed_orf_list.fasta_uniq.fa < Current_species_T2R.txt > Functionnal_T2Rs_multifasta_singleexon.fa

#Remove sequences with ambigous nucleotides
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_T2Rs_multifasta_singleexon.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Functionnal_T2Rs_multifasta_singleexon.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_T2Rs_multifasta_singleexon.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > unclear_Functionnal_T2Rs_multifasta_singleexon.fa



#Translate coding sequences
transeq clear_Functionnal_T2Rs_multifasta_singleexon.fa clear_Functionnal_T2Rs_multifasta_singleexon.prot ; sed -i 's/_1$//g' clear_Functionnal_T2Rs_multifasta_singleexon.prot



###Extract the coordinates of functionnals ORs found and perform a tblastn with found genes

grep ">" Functionnal_T2Rs_multifasta_singleexon.fa | sed 's/>//g' | sed 's/-/	/g' > Coordinates_Functionnal_T2R.txt
if [ `wc -l < Coordinates_Functionnal_T2R.txt` -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_Functionnal_T2R.txt ; fi


tblastn -query clear_Functionnal_T2Rs_multifasta_singleexon.prot -db $genome -evalue $evalue -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out tblastn_functionnal_t2r_vs_genome.tblastn -num_threads $number_of_thread

#filter the tblastn to remove problematic sequences
cut -f1 tblastn_functionnal_t2r_vs_genome.tblastn  | sort | uniq -c | sed 's/^ *//g'| sed 's/ /,/g' | awk -F, ' ($1 < 2000) ' | cut -f2 -d "," > good_tblastn_id ; for i in `cat good_tblastn_id` ; do grep "$i" tblastn_functionnal_t2r_vs_genome.tblastn >> Filtered_tblastn_functionnal_t2r_vs_genome.tblastn ; done


cat T2R_vs_Genome.blastn Filtered_tblastn_functionnal_t2r_vs_genome.tblastn > tblastn_functionnal_and_known_t2r_vs_genome.tblastn


cat $T2R_database clear_Functionnal_T2Rs_multifasta_singleexon.prot > Complete_T2R_db.prot



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Extraction of singe exon T2R pseudogenes  ####################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


#Launch a Rscript. This will find regions corresponding to single exon 
#pseudogenes/ truncated genes by parsing the tblastn results and removing already analyzed regions

Rscript $scripts_location/R_script_numero2_2.R #This Rscript generate "Pseudo_truncated_coordinates.tsv" with 4 columns : scaffold, coord_start, coord_end, query, non-extanded start, non-extanded stop

#extract regions containing putative T2R pseudogenes

cut -f1 Pseudo_truncated_coordinates.tsv > scaffolds.txt
cut -f2,3 Pseudo_truncated_coordinates.tsv | sed 's/	/-/g' > coordinates.txt
paste -d ":" scaffolds.txt coordinates.txt > Pseudogenes_regions.tsv

xargs samtools faidx $genome < Pseudogenes_regions.tsv > Pseudogenes_regions.fa


#Retain only besthit that have a T2R gene in atleast the 3 first hits
blastx -query Pseudogenes_regions.fa -db $blast_database -max_target_seqs 3 -outfmt '6 qseqid sseqid' -out blastx_blast_regions.tsv -num_threads $number_of_thread
grep "T2R-Receptor\|Taste-receptor" blastx_blast_regions.tsv | cut -f1 | sort | uniq > T2R_best_hits_regions.tsv
xargs samtools faidx $genome < T2R_best_hits_regions.tsv > T2R_best_hits_regions.fa


IFS=$'\n'
for line in `cat T2R_best_hits_regions.tsv` ; do 
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
	samtools faidx Complete_T2R_db.prot $query > Gene.prot

	blastx -query Genomic_region.fa -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastx_result -max_target_seqs 1 -num_threads 10
	 
	if grep -q -i "T2R-Receptor\|Taste-receptor" blastx_result ; then

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




#IFS=$'\n' ; for line in `cat Pseudo_truncated_coordinates.tsv` ; do scaffold=`echo $line | cut -f1` ; start=`echo $line | cut -f2` ; stop=`echo $line | cut -f3` ; query=`echo $line | cut -f4` ; scaffold_length=`grep "^$scaffold " $genome.fai | cut -f2` ; true_start=`echo $line | cut -f5` ; true_end=`echo $line | cut -f6` ; strand=`echo $line | cut -f7` ; stop_codon_state="FALSE" ; edge_state="FALSE" ; frameshift_state="FALSE" ; samtools faidx $genome $scaffold:$start-$stop > Genomic_region.fa ; makeblastdb -in Genomic_region.fa -dbtype nucl ; samtools faidx Complete_T2R_db.prot $query > Gene.prot ; blastx -query Genomic_region.fa -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastx_result -max_target_seqs 1 -num_threads 10 ; if grep -q -i "T2R-Receptor\|Taste-receptor" blastx_result ; then tblastn -query Gene.prot -db Genomic_region.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe sseq' > tblast_result_stop.tsv ; tblastn -gapopen 32767 -gapextend 32767 -query Gene.prot -db Genomic_region.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe sseq' > tblast_result_fs.tsv ; cut -f12 tblast_result_fs.tsv > Frame_detection.txt ; cut -f13 tblast_result_stop.tsv > Stop_detection.txt ; number_hsp=`cat Frame_detection.txt | wc -l`  ; number_diff_frames=`cat Frame_detection.txt | sort | uniq | wc -l` ; if [ "$number_diff_frames" -gt '1' ] ; then frameshift_state="TRUE" ; fi ; number_stops=`grep -o "\*" Stop_detection.txt | wc -l` ; if [ "$number_stops" -ge '1' ] ; then stop_codon_state="TRUE" ; fi ; if [ "$true_start" -le '100' ] ; then edge_state="TRUE" ; fi  ; diff_lengths=$((scaffold_length - true_end)) ; if [ "$diff_lengths" -le '100' ] ; then edge_state="TRUE" ; fi ; extanded_start_coord=$((start - 100)) ; extanded_end_coord=$((stop + 100)) ; consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1` ; if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi ; if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi  ; if [ $strand == "+" ] ; then samtools faidx $genome $scaffold:$true_start-$true_end > temporary_rslt.fa ; grep -v ">" temporary_rslt.fa > temporary_rslt.txt ; header_name=`echo "$scaffold-$true_start-$true_end---1_exon-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'` ; sed -e "1i>$header_name\\"  temporary_rslt.txt > temporary_rslt_renamed.fa ; sed '/^[[:space:]]*$/d' temporary_rslt_renamed.fa >> Pseudogenes_single_exon.fa ; $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 Gene.prot Genomic_region.fa > Exonerate_single_exon_pseudo ; awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' Exonerate_single_exon_pseudo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa ; transeq predicted_cds.fa predicted_cds.prot ; rm predicted_cds.fa ; sed "s/>.*/>$header_name/g" predicted_cds.prot >> Frameshiftless_Pseudogenes_single_exon.prot ; else samtools faidx $genome $scaffold:$true_start-$true_end > temporary_rslt.fa ; revseq temporary_rslt.fa temporary_rslt_rev.fa ; grep -v ">" temporary_rslt_rev.fa > temporary_rslt_rev.txt ; header_name=`echo "$scaffold-$true_start-$true_end---1_exon-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g'` ; sed -e "1i>$header_name\\" temporary_rslt_rev.txt > temporary_rslt_renamed.fa ; sed '/^[[:space:]]*$/d' temporary_rslt_renamed.fa >> Pseudogenes_single_exon.fa ; $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 Gene.prot Genomic_region.fa > Exonerate_single_exon_pseudo  ; awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' Exonerate_single_exon_pseudo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa ; transeq predicted_cds.fa predicted_cds.prot ; rm predicted_cds.fa ; sed "s/>.*/>$header_name/g" predicted_cds.prot >> Frameshiftless_Pseudogenes_single_exon.prot ; fi ;fi ; done

# Remove pseudogenes with ambigous nucleotides or a length less than 150nt
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_single_exon.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > Pseudogenes_single_exon_clean.fa 
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_single_exon.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > Pseudogenes_single_exon_unclean.fa 
grep ">" Pseudogenes_single_exon_clean.fa  | sed 's/>//g' > clean_id ; for id in `cat clean_id` ; do if grep -q "$id" Frameshiftless_Pseudogenes_single_exon.prot ; then echo "ok" ; else sed -i "s/$id//g" clean_id ; fi ; done ; xargs samtools faidx Frameshiftless_Pseudogenes_single_exon.prot < clean_id > Frameshiftless_Pseudogenes_single_exon_clean.prot
grep ">" Pseudogenes_single_exon_unclean.fa  | sed 's/>//g' > unclean_id ; for id in `cat unclean_id` ; do if grep -q "$id" Frameshiftless_Pseudogenes_single_exon.prot ; then echo "ok" ; else sed -i "s/$id//g" unclean_id ; fi ; done ; xargs samtools faidx Frameshiftless_Pseudogenes_single_exon.prot < unclean_id > Frameshiftless_Pseudogenes_single_exon_unclean.prot


seqkit seq -m 150 Pseudogenes_single_exon_clean.fa  > Pseudogenes_single_exon_clean_s200.fa 


###NC_030681.2-160010718-160011479

#Verify that these gene fragments correspong to T2R genes with a blastx
blastx -query Pseudogenes_single_exon_clean_s200.fa -db $blast_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastx_result -max_target_seqs 1 -num_threads $number_of_thread
grep -i "T2R-Receptor\|Taste-receptor" blastx_result | cut -f1 | sort | uniq > single_exon_pseudos_from_blast.list

fasta_formatter -i Pseudogenes_single_exon_clean_s200.fa > reformat_Pseudogenes_single_exon_clean_s200.fa ; mv reformat_Pseudogenes_single_exon_clean_s200.fa Pseudogenes_single_exon_clean_s200.fa
xargs samtools faidx Pseudogenes_single_exon_clean_s200.fa < single_exon_pseudos_from_blast.list > Verified_Pseudogenes_single_exon_clean_s200.fa

for id in `cat single_exon_pseudos_from_blast.list` ; do if grep -q "$id" Frameshiftless_Pseudogenes_single_exon_clean.prot ; then echo "ok" ; else sed -i "s/$id//g" single_exon_pseudos_from_blast.list ; fi ; done
xargs samtools faidx Frameshiftless_Pseudogenes_single_exon_clean.prot < single_exon_pseudos_from_blast.list > Verified_Frameshiftless_Pseudogenes_single_exon_clean.fa


#Verify that pseudogenes correspond to T2R genes with a ML tree


if [ `grep -c ">" Verified_Frameshiftless_Pseudogenes_single_exon_clean.fa` -ge 1 ] ; then 
	mafft --add Verified_Frameshiftless_Pseudogenes_single_exon_clean.fa --keeplength $scripts_location/T2R_plus_outgroup.aln > Putative_t2r_plus_known_t2r_plus_outgroup.prot.aln
	iqtree -s Putative_t2r_plus_known_t2r_plus_outgroup.prot.aln -st AA -nt $number_of_thread -m JTT+F+I+G4 -redo -n 200
	Rscript $scripts_location/Tree_parser.R 
	xargs samtools faidx Verified_Pseudogenes_single_exon_clean_s200.fa < Current_species_T2R.txt > Final_Pseudogenes.fa
else 
	echo "" > Final_Pseudogenes.fa
fi



#Rename genes

grep ">" clear_Functionnal_T2Rs_multifasta_singleexon.fa | sed 's/>//g' > single_exon_id.txt
sed -e 's/$/---1_exons/' single_exon_id.txt > single_exon_id_edit.txt
paste -d "\t" single_exon_id.txt single_exon_id_edit.txt > renaming_single_exon.txt
perl $scripts_location/rename_fasta.pl renaming_single_exon.txt clear_Functionnal_T2Rs_multifasta_singleexon.fa > temporary.fasta ; mv temporary.fasta clear_Functionnal_T2Rs_multifasta_singleexon.fa

cp clear_Functionnal_T2Rs_multifasta_singleexon.fa Combined_Functionnal_T2Rs.fa


#Check if complete genes have 7tm domain determined with phobius or TMHMM


#First check with phobius
transeq Combined_Functionnal_T2Rs.fa Combined_Functionnal_T2Rs.prot ; sed -i 's/_1$//g' Combined_Functionnal_T2Rs.prot #translate CDS
perl $scripts_location/phobius/phobius.pl -long Combined_Functionnal_T2Rs.prot > Phobius_verification.txt #run phonius in long mode
grep ">" Combined_Functionnal_T2Rs.prot | sed 's/>//g' > gene_id.txt #extract cds id

for gene in `cat gene_id.txt` ; do nb_transm=`sed '/'"$gene"'/,/\/\//!d;/\/\//q' Phobius_verification.txt | grep "TRANSMEM" | wc -l` ; echo "$gene,$nb_transm" ; done > Gene_NbTm.tsv
awk 'BEGIN{FS=",";OFS=","}($2>=7){print $1;}' Gene_NbTm.tsv > Phobius_genes_with_7tm.txt
awk 'BEGIN{FS=",";OFS=","}($2<7){print $1;}' Gene_NbTm.tsv > Phobius_genes_without_7tm.txt

#Now, with TMHMM

$scripts_location/tmhmm-2.0c/bin/tmhmm Combined_Functionnal_T2Rs.prot > tmhmm_verification.txt
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

perl $scripts_location/rename_fasta.pl renaming_file_tm.txt Combined_Functionnal_T2Rs.fa > temporary.fasta ; mv temporary.fasta Combined_Functionnal_T2Rs.fa


#merge the two file containing ambigous sequences

cat Pseudogenes_single_exon_unclean.fa unclear_Functionnal_T2Rs_multifasta_singleexon.fa > Ambigous_T2R.fasta


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

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Ambigous_T2R.fasta | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Ambigous_T2R_uniq.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Final_Pseudogenes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Final_Pseudogenes_uniq.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Combined_Functionnal_T2Rs.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Combined_Functionnal_T2Rs_uniq.fa



nb_seq=`grep -c ">" Ambigous_T2R_uniq.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" Ambigous_T2R_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_ambigous_final.tsv
fi

nb_seq=`grep -c ">" Final_Pseudogenes_uniq.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" Final_Pseudogenes_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_pseudogenes_final.tsv
fi

nb_seq=`grep -c ">" Combined_Functionnal_T2Rs_uniq.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" Combined_Functionnal_T2Rs_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_genes_final.tsv
fi




Rscript $scripts_location/Remove_redundancy.R



IFS=$'\n'


for line in `cat best_genes_functionnal.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Combined_Functionnal_T2Rs_uniq.fa | sed 's/>//g' >> functionnal_to_keep.txt ; done
for line in `cat best_genes_ambigous.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Ambigous_T2R_uniq.fa | sed 's/>//g' >> ambigous_to_keep.txt ; done
for line in `cat best_genes_pseudogenes.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" Final_Pseudogenes_uniq.fa | sed 's/>//g' >> pseudogenes_to_keep.txt ; done


xargs samtools faidx Combined_Functionnal_T2Rs_uniq.fa < functionnal_to_keep.txt > FINAL_Functionnal_T2R.fa
xargs samtools faidx Ambigous_T2R_uniq.fa < ambigous_to_keep.txt > FINAL_Ambigous_T2R.fa
xargs samtools faidx Final_Pseudogenes_uniq.fa < pseudogenes_to_keep.txt > FINAL_Pseudogenes_T2R.fa




#Also classify genes without 7tm as pseudogenes (if the option is TRUE)


if [ $tm_filter == "TRUE" ] ; then

	grep ">" FINAL_Functionnal_T2R.fa | grep "phobius\|tmhmm" | sed 's/>//g' > 7tm_genes
	grep ">" FINAL_Functionnal_T2R.fa | grep -v "phobius\|tmhmm" | sed 's/>//g' > non_7tm_genes
	xargs samtools faidx FINAL_Functionnal_T2R.fa < 7tm_genes > FINAL_Functionnal_T2R_7tm.fa 
	xargs samtools faidx FINAL_Functionnal_T2R.fa < non_7tm_genes >> FINAL_Pseudogenes_T2R.fa

else
	cp FINAL_Functionnal_T2R.fa FINAL_Functionnal_T2R_7tm.fa 
fi



## FINAL FILTER

transeq FINAL_Functionnal_T2R_7tm.fa FINAL_Functionnal_T2R_7tm.prot ; sed -i 's/_1$//g' FINAL_Functionnal_T2R_7tm.prot
blastp -query FINAL_Functionnal_T2R_7tm.prot -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastp_result_custom -max_target_seqs 1 -num_threads $number_of_thread
grep "T2R\|Taste-receptor" blastp_result_custom | cut -f1 | sort | uniq  > good_sequences
xargs samtools faidx FINAL_Functionnal_T2R_7tm.fa < good_sequences > temp.fa
mv temp.fa FINAL_Functionnal_T2R_7tm.fa ; rm good_sequences ; rm *.fai


blastx -query FINAL_Pseudogenes_T2R.fa -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastx_result_custom -max_target_seqs 1 -num_threads $number_of_thread
grep "T2R\|Taste-receptor" blastx_result_custom | cut -f1 | sort | uniq  > good_sequences_p
xargs samtools faidx FINAL_Pseudogenes_T2R.fa < good_sequences_p > temp_p.fa
mv temp_p.fa FINAL_Pseudogenes_T2R.fa ; rm good_sequences_p ; rm *.fai 


blastx -query FINAL_Ambigous_T2R.fa -db $blast_database -outfmt '6 qseqid sseqid evalue' -out blastx_result_custom -max_target_seqs 1 -num_threads $number_of_thread
grep "T2R\|Taste-receptor" blastx_result_custom | cut -f1 | sort | uniq  > good_sequences_p
xargs samtools faidx FINAL_Ambigous_T2R.fa < good_sequences_p > temp_p.fa
mv temp_p.fa FINAL_Ambigous_T2R.fa ; rm good_sequences_p ; rm *.fai



#We have three final file :
#All potentially functionnal genes : FINAL_Functionnal_T2R_7tm.fa
#Probably pseudogenes or edge genes : FINAL_Pseudogenes_T2R.fa
#Ambigous sequences : FINAL_Ambigous_T2R.fa

nb_functionnal=`grep -c ">" FINAL_Functionnal_T2R_7tm.fa`
nb_pseudo_edge=`grep -c ">" FINAL_Pseudogenes_T2R.fa`
nb_ambigous=`grep -c ">" FINAL_Ambigous_T2R.fa`

echo "Search of T2R is finished. There are $nb_functionnal potentially functionnal genes, $nb_pseudo_edge pseudogenes or fragments and $nb_ambigous ambigous sequences"

echo "$nb_functionnal	$nb_pseudo_edge	$nb_ambigous" > Results_NbF_NbP_NbA_summary.txt



dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"

