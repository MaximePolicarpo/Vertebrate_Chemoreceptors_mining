#!/bin/bash


#Put the genome as first argument

genome=$1


#Prepare the genome and perform tblastn of known OR genes against it

sed -i "s/ .*//g" $genome ; sed -i "s/a/A/g" $genome ; sed -i "s/t/T/g" $genome ; sed -i "s/g/G/g" $genome ; sed -i "s/c/C/g" $genome
makeblastdb -in $genome -dbtype nucl
samtools faidx $genome

tblastn -query /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot.0 -db $genome -evalue 1e-5 -outfmt 6 -out OR_vs_Genome_0.blastn -num_threads 10
tblastn -query /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot.1 -db $genome -evalue 1e-5 -outfmt 6 -out OR_vs_Genome_1.blastn -num_threads 10
tblastn -query /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot.2 -db $genome -evalue 1e-5 -outfmt 6 -out OR_vs_Genome_2.blastn -num_threads 10
tblastn -query /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot.3 -db $genome -evalue 1e-5 -outfmt 6 -out OR_vs_Genome_3.blastn -num_threads 10
tblastn -query /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot.4 -db $genome -evalue 1e-5 -outfmt 6 -out OR_vs_Genome_4.blastn -num_threads 10


cat OR_vs_Genome_0.blastn OR_vs_Genome_1.blastn OR_vs_Genome_2.blastn OR_vs_Genome_3.blastn OR_vs_Genome_4.blastn > OR_vs_Genome.blastn



#Now lets perform a tblastn with a more stringent evalue to find truncated, pseudogenes and edge genes and maybe multiexon OR genes
tblastn -query /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot.0 -db $genome -evalue 1e-20 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out tblastn_functionnal_or_vs_genome_0.blastn -num_threads 10
tblastn -query /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot.1 -db $genome -evalue 1e-20 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out tblastn_functionnal_or_vs_genome_1.blastn -num_threads 10
tblastn -query /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot.2 -db $genome -evalue 1e-20 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out tblastn_functionnal_or_vs_genome_2.blastn -num_threads 10
tblastn -query /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot.3 -db $genome -evalue 1e-20 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out tblastn_functionnal_or_vs_genome_3.blastn -num_threads 10
tblastn -query /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot.4 -db $genome -evalue 1e-20 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out tblastn_functionnal_or_vs_genome_4.blastn -num_threads 10

cat tblastn_functionnal_or_vs_genome_0.blastn tblastn_functionnal_or_vs_genome_1.blastn tblastn_functionnal_or_vs_genome_2.blastn tblastn_functionnal_or_vs_genome_3.blastn tblastn_functionnal_or_vs_genome_4.blastn > tblastn_functionnal_or_vs_genome.tblastn


#Extract non-overlapping best-hit sequences (first round of tblastn with the evalue of 1e-10)

/home/policarpo/miniconda3/bin/Rscript /mnt/35To/policarpo/2021_OR_Finder/Database/R_script_numero1.R 


#For each best-hits we will extract ORFs that are longer than 750 bp

rm Putative_or_multifasta_raw.fa

IFS=$'\n'

for line in `cat Best_hits_filtered.tsv` ; do
	
	echo "$line" > test


	scaffold=`cut -f1 test`
	initial_start=`cut -f2 test` 
	initial_stop=`cut -f3 test` 

	new_3prime=$((initial_start - 1000)) #extend the best-hit 1000bp upstream an downtream to make sure we get the ATG and TAG
	new_5prime=$((initial_stop + 1000)) 

	if [ "$new_3prime" -lt '1' ] ; then  #If the extension go below 1 then we put it to 1 (scaffold beginning)
		new_3prime=1
	fi


	samtools faidx $genome $scaffold:$new_3prime-$new_5prime > sequence.fa 
	getorf -sequence sequence.fa -outseq orf_list -minsize 750 -find 3 
	sed -i "s/>/>$scaffold:/g" orf_list ; sed -i 's/ //g' orf_list
	cat orf_list >> Putative_or_multifasta_raw.fa

done

sed -i 's/(REVERSESENSE)/_reverse/g' Putative_or_multifasta_raw.fa
grep ">" Putative_or_multifasta_raw.fa | sed 's/>//g' > oldfastaheaders

sed 's/:/	/g' oldfastaheaders | sed 's/-/	/g' | sed 's/_[0-9]\[/	/g' | sed 's/\]//g' | sed 's/_reverse/	reverse/g' > temporary


#Lets rename sequences so that it has the format >ScaffoldName-StartCoor-EndCoor

IFS=$'\n'

for line in `cat temporary`
	do
		echo "$line" > temp 

		scaffold=`cut -f1 temp`

		if grep -q "reverse" temp ; then
			coord_start=`awk 'BEGIN {FS="	"} {sum+=$2+$5-1} END {print sum}' temp`
			coord_end=`awk 'BEGIN {FS="	"} {sum+=$2+$4-1} END {print sum}' temp`
		else
			coord_start=`awk 'BEGIN {FS="	"} {sum+=$2+$4-1} END {print sum}' temp`
			coord_end=`awk 'BEGIN {FS="	"} {sum+=$2+$5-1} END {print sum}' temp`

		fi


		echo "$scaffold-$coord_start-$coord_end" >> newfastaheaders
done


paste -d "	" oldfastaheaders newfastaheaders > renaming_file

perl /mnt/35To/policarpo/2021_OR_Finder/Database/rename_fasta.pl renaming_file Putative_or_multifasta_raw.fa > Putative_or_multifasta.fa

#Due to the -1000/+1000 extension, some orf could be found twice so we remove identical sequences

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Putative_or_multifasta.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Putative_or_multifasta_uniq.fa

#Lets translate DNA to Prot

transeq Putative_or_multifasta_uniq.fa Putative_or_multifasta.prot  ; sed -i 's/_1$//g' Putative_or_multifasta.prot 

#Lets blastp our sequences against uniprot database to filter real OR genes from artifacts

blastp -query Putative_or_multifasta.prot -db /mnt/35To/policarpo/DataBase_Uniprot/uniprot_sprot.fasta -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 20
grep -i "olfactory\|odorant" blastp_result | cut -f1 | sort | uniq > OR_from_blast.list
for i in `cat OR_from_blast.list` ; do samtools faidx Putative_or_multifasta.prot $i >> Fasta_OR_from_blast.prot ; done 



#Lets align our sequences and perform a ML tree with known OR genes
mafft --add Fasta_OR_from_blast.prot --keeplength /mnt/35To/policarpo/2021_OR_Finder/Database/Phylogenetic_tree_Database/Template_alignment > Putative_or_plus_known_or_plus_outgroup.prot.aln
/home/policarpo/iqtree-1.6.12-Linux/bin/iqtree -s Putative_or_plus_known_or_plus_outgroup.prot.aln -st AA -nt 50 -bb 1000 -m JTT+F+I+G4 -redo



#We use a script to remove outgroup sequences that are not OR
cp /mnt/35To/policarpo/2021_OR_Finder/Database/Known_OR_genes_TYPEI.id ./
cp /mnt/35To/policarpo/2021_OR_Finder/Database/Known_OR_genes_TYPEII.id ./
/home/policarpo/miniconda3/bin/Rscript /mnt/35To/policarpo/2021_OR_Finder/Database/Tree_parser.R


#Remaining sequences are OR. Lets remove 100% identical sequences and remove complete sequences but with ambigous nucleotides
for i in `cat Current_species_OR.txt`  ; do samtools faidx Putative_or_multifasta_uniq.fa $i ; done > Functionnal_ORs_multifasta.fa
/home/policarpo/cdhit-master/cd-hit -i Functionnal_ORs_multifasta.fa -o Functionnal_ORs_multifasta_cdhit.fa -c 1
python2 /mnt/35To/policarpo/Projet_Confinement_OR/clearer_ambigous_nt.py Functionnal_ORs_multifasta_cdhit.fa 0 0 


transeq clear_Functionnal_ORs_multifasta_cdhit.fa clear_Functionnal_ORs_multifasta_cdhit.prot ; sed -i 's/_1$//g' clear_Functionnal_ORs_multifasta_cdhit.prot


IFS=$'\n'

grep ">" clear_Functionnal_ORs_multifasta_cdhit.fa | sed 's/>//g' > clears_ID.list
grep ">" Functionnal_ORs_multifasta_cdhit.fa | sed 's/>//g' > unclear_IDs.list

for i in `cat unclear_IDs.list` 
	do 
		if grep -q "$i" clears_ID.list

			then echo "clear"
		
		else 
			echo "$i" >> Ambiguous_sequences.list

		fi
done


for i in `cat Ambiguous_sequences.list` ; do samtools faidx Functionnal_ORs_multifasta_cdhit.fa $i >> Ambigous_sequences_multifasta.fa ; done




#extract the coordinates of functionnals ORs found. Then use the script to perform the tblastn with an evalue of 1e-20

grep ">" Functionnal_ORs_multifasta.fa | sed 's/>//g' | sed 's/-/	/g' > Coordinates_Functionnal_ORS.txt

tblastn -query clear_Functionnal_ORs_multifasta_cdhit.prot -db $genome -evalue 1e-20 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out tblastn_functionnal_or_vs_genome_current_species.tblastn -num_threads 10

cat tblastn_functionnal_or_vs_genome_0.blastn tblastn_functionnal_or_vs_genome_1.blastn tblastn_functionnal_or_vs_genome_2.blastn tblastn_functionnal_or_vs_genome_3.blastn tblastn_functionnal_or_vs_genome_4.blastn tblastn_functionnal_or_vs_genome_current_species.tblastn > tblastn_functionnal_or_vs_genome.tblastn


/home/policarpo/miniconda3/bin/Rscript /mnt/35To/policarpo/2021_OR_Finder/Database/R_sript_merge_multiexon.R


### Search for potential multiple exon regions



cat /mnt/35To/policarpo/2021_OR_Finder/Database/Database_Vertebrate_OR_cdhit_80.prot clear_Functionnal_ORs_multifasta_cdhit.prot > ALL_functional_OR.prot


IFS=$'\n'

rm exonerate.test
for line in `cat Potential_multiple_exon_regions.tsv` ; do 
	echo "$line" > Current_line.txt

	scaffold=`cut -f1 Current_line.txt`
	begin=`cut -f2 Current_line.txt`
	end=`cut -f3 Current_line.txt`
	query=`cut -f4 Current_line.txt`

	samtools faidx $genome $scaffold:$begin-$end > scaffold.fa

	samtools faidx ALL_functional_OR.prot $query > protein.prot

	/home/casane/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome --minintron 80 --maxintron 4000 protein.prot scaffold.fa >> exonerate.test

done


grep "vulgar:" exonerate.test  > vulgar_lines.txt



IFS=$'\n'

rm Pseudo_multi_exon.fa
rm Multiple_exon_functional_genes.fa
rm Truncated_multi_exon.fa
rm Edge_multi_exon.fa

for line in `cat vulgar_lines.txt` ; do 

	rm orf_list.fai
	rm Current_sequence.fa.fai
	rm Current_sequence.fa


	echo "$line" > current_line_vulgar.txt
	sed 's/vulgar: //g' current_line_vulgar.txt | cut -f2- -d " " > current_line_vulgar_cut.txt


	#Extract query infos

	query=`cut -f2 -d " " current_line_vulgar.txt`
	query_start=`cut -f1 -d " " current_line_vulgar_cut.txt`
	query_end=`cut -f2 -d " " current_line_vulgar_cut.txt`


	#Extract scaffold infos

	scaffold=`cut -f4 -d " " current_line_vulgar_cut.txt`
	scaffold_name=`echo "$scaffold" | sed 's/:/ /g' | sed 's/-/ /g' | cut -f1 -d " "`
	scaffold_start=`echo "$scaffold" | sed 's/:/ /g' | sed 's/-/ /g' | cut -f2 -d " "`
	scaffold_end=`echo "$scaffold" | sed 's/:/ /g' | sed 's/-/ /g' | cut -f3 -d " "`

	scaffold_length=`grep "^$scaffold_name	" $genome.fai | cut -f2`
	contig_end=$((scaffold_length-30))

	#Extract target infos

	target_start=`cut -f5 -d " " current_line_vulgar_cut.txt`
	target_end=`cut -f6 -d " " current_line_vulgar_cut.txt`
	echo "$target_start" > coordinates.txt
	echo "$target_end" >> coordinates.txt
	target_start=`sort -n coordinates.txt | head -1`
	target_end=`sort -n coordinates.txt | tail -1`

	strand=`cut -f7 -d " " current_line_vulgar_cut.txt`

	#Extand a bit the target sequence

	new_start=$((target_start - 300))
	new_end=$((target_end + 300))


	#transform it on true genome coordinate

	genomic_coord_start=$((scaffold_start + new_start))

	if [ $genomic_coord_start -lt 2 ] ; then genomic_coord_start=1 ; fi

	genomic_coord_end=$((scaffold_start + new_end))

	query_length_aln=$((query_end - query_start))

	line_current=`cat current_line_vulgar.txt`

	if [ $query_length_aln -gt 200 ] ; then

		if grep -q " F " current_line_vulgar_cut.txt ; then

			samtools faidx ALL_functional_OR.prot $query > current_protein.prot 
			samtools faidx $genome $scaffold_name:$genomic_coord_start-$genomic_coord_end > Current_sequence.fa
			sed -i 's/:/-/g' Current_sequence.fa

			/home/casane/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome --minintron 80 --maxintron 4000 --ryo "%tcs" --showtargetgff True --bestn 1 current_protein.prot Current_sequence.fa > current_exonerate.exo


			exon_number=`awk '/##gff-version 2/ {p=1}; p; /# --- END OF GFF DUMP ---/ {p=0}' current_exonerate.exo | cut -f3 | grep -v "#" | grep -c "exon"`


			if [ "$exon_number" -gt '1' ] ; then

				awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' current_exonerate.exo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa

				sed -i "s/>.*/>$scaffold_name:$genomic_coord_start-$genomic_coord_end-$exon_number-exons/g" predicted_cds.fa

				cat predicted_cds.fa >> Pseudo_multi_exon.fa
			fi



		else

			samtools faidx ALL_functional_OR.prot $query > current_protein.prot 
			samtools faidx $genome $scaffold_name:$genomic_coord_start-$genomic_coord_end > Current_sequence.fa
			sed -i 's/:/-/g' Current_sequence.fa
			/home/casane/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome --minintron 80 --maxintron 4000 --ryo "%tcs" --showtargetgff True --bestn 1 current_protein.prot Current_sequence.fa > current_exonerate.exo

			target_range_1=`grep -m1 "Target range:" current_exonerate.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
			target_range_2=`grep -m1 "Target range:" current_exonerate.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`

			echo "$target_range_1" > Target_range.txt
			echo "$target_range_2" >> Target_range.txt

			target_range_start=`sort -n Target_range.txt | head -1`
			target_range_end=`sort -n Target_range.txt | tail -1`

			exon_number=`awk '/##gff-version 2/ {p=1}; p; /# --- END OF GFF DUMP ---/ {p=0}' current_exonerate.exo | cut -f3 | grep -v "#" | grep -c "exon"`


			if [ $strand == "+" ] && [ "$exon_number" -gt '1' ] ; then 



				seqname=`grep ">" Current_sequence.fa | sed 's/>//'`

				samtools faidx Current_sequence.fa $seqname:1-$target_range_start > Extend_three_prime.fa
				grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt



				read_target_end=$((target_range_end + 1))
				read_extansion_end=$((target_range_end + 301))

				samtools faidx Current_sequence.fa $seqname:$read_target_end-$read_extansion_end > Extend_five_prime.fa
				grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt


				awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' current_exonerate.exo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
				grep -v ">" predicted_cds.fa > predicted_cds.txt

				cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_seq.fa


				sed -i '1 i\>Complete_seq' Complete_extanded_seq.fa

				getorf -sequence Complete_extanded_seq.fa -outseq orf_list -minsize 750 -find 3

				if [ `grep -c ">" orf_list` -ge 1 ] ; then

					
					sed -i 's/ //g' orf_list
					sequence_name=`grep -m1 ">" orf_list | sed 's/>//g'`

					samtools faidx orf_list $sequence_name > myseq.fa
					transeq myseq.fa myseq.prot

					/home/casane/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --minintron 80 --maxintron 4000 --ryo "%tcs" --showtargetgff True --bestn 1 myseq.prot Current_sequence.fa > current_exonerate_parsed.exo


					target_range_final_1=`grep -m1 "Target range:" current_exonerate_parsed.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					target_range_final_2=`grep -m1 "Target range:" current_exonerate_parsed.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
					
					coord_start_gene=$((genomic_coord_start + target_range_final_1))
					coord_end_gene=$((genomic_coord_start + target_range_final_2))

					sed "s/>.*/>$scaffold_name:$coord_start_gene-$coord_end_gene-$exon_number-exons/g" myseq.fa >> Multiple_exon_functional_genes.fa

				elif [ "$exon_number" -gt '1' ] ; then

					awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' current_exonerate.exo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
					sed -i "s/>.*/>$scaffold_name:$genomic_coord_start-$genomic_coord_end-$exon_number-exons/g" predicted_cds.fa

					transeq predicted_cds.fa predicted_cds.prot ; sed -i 's/_1$//g' predicted_cds.prot

					if grep -q "\*" predicted_cds.prot ; then 
						cat predicted_cds.fa >> Pseudo_multi_exon.fa
					elif [ $genomic_coord_start -le '30' ] || [ $genomic_coord_end -gt $contig_end ] ; then 
						cat predicted_cds.fa >> Edge_multi_exon.fa
					else
						cat predicted_cds.fa >> Truncated_multi_exon.fa

					fi

				fi

			elif [ $strand == "-" ] && [ "$exon_number" -gt '1' ] ; then

				seqname=`grep ">" Current_sequence.fa | sed 's/>//'`

				read_target_end=$((target_range_end + 1))
				read_extansion_end=$((target_range_end + 301))

				samtools faidx Current_sequence.fa $seqname:$read_target_end-$read_extansion_end > Extend_three_prime.fa
				grep -v ">" Extend_three_prime.fa > Extend_three_prime.txt

				samtools faidx Current_sequence.fa $seqname:1-$target_range_start > Extend_five_prime.fa
				grep -v ">" Extend_five_prime.fa > Extend_five_prime.txt

				awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' current_exonerate.exo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
				revseq predicted_cds.fa predicted_cds_rev.fa
				grep -v ">" predicted_cds_rev.fa > predicted_cds_rev.txt

				cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_seq.fa

				sed -i '1 i\>Complete_seq' Complete_extanded_seq.fa

				getorf -sequence Complete_extanded_seq.fa -outseq orf_list -minsize 750 -find 3

				if [ `grep -c ">" orf_list` -ge 1 ] ; then

					sed -i 's/ //g' orf_list
					sequence_name=`grep -m1 ">" orf_list | sed 's/>//g'`

					samtools faidx orf_list $sequence_name > myseq.fa
					transeq myseq.fa myseq.prot
				

					/home/casane/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --minintron 50 --maxintron 4000 --ryo "%tcs" --showtargetgff True --bestn 1 myseq.prot Current_sequence.fa > current_exonerate_parsed.exo


					target_range_final_1=`grep -m1 "Target range:" current_exonerate_parsed.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
					target_range_final_2=`grep -m1 "Target range:" current_exonerate_parsed.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`

					coord_start_gene=$((genomic_coord_start + target_range_final_2))
					coord_end_gene=$((genomic_coord_start + target_range_final_1))

					sed "s/>.*/>$scaffold_name:$coord_start_gene-$coord_end_gene-$exon_number-exons/g" myseq.fa >> Multiple_exon_functional_genes.fa

				elif [ "$exon_number" -gt '1' ] ; then

					awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' current_exonerate.exo | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
					sed -i "s/>.*/>$scaffold_name:$genomic_coord_start-$genomic_coord_end-$exon_number-exons/g" predicted_cds.fa

					if grep -q "\*" predicted_cds.prot ; then 
						cat predicted_cds.fa >> Pseudo_multi_exon.fa
					elif
						[ $genomic_coord_start -le '30' ] || [ $genomic_coord_end -gt $contig_end ] ; then
							cat predicted_cds.fa >> Edge_multi_exon.fa
					else
						cat predicted_cds.fa >> Truncated_multi_exon.fa

					fi

				fi

			fi

		fi

	fi


done



grep ">" Multiple_exon_functional_genes.fa | sed 's/>//g' | sed 's/:/	/g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_Functionnal_ORS.txt
grep ">" Pseudo_multi_exon.fa | sed 's/>//g' | sed 's/:/	/g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_Functionnal_ORS.txt
grep ">" Truncated_multi_exon.fa | sed 's/>//g' | sed 's/:/	/g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_Functionnal_ORS.txt
grep ">" Edge_multi_exon.fa | sed 's/>//g' | sed 's/:/	/g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_Functionnal_ORS.txt



#Lets do a filter to only retain regions correspond to edge, truncated or pseudogenes



/home/policarpo/miniconda3/bin/Rscript /mnt/35To/policarpo/2021_OR_Finder/Database/R_script_numero2_2.R



#Now we are goind to make a differenciation between pseudogenes and truncated genes. pseudogenes have one LoF mutation in the tblastn result #
#Truncated genes are genes that are <30bp close to the beginning or end of a contig with no LoF

IFS=$'\n'

for i in `cat Pseudo_truncated_coordinates.tsv` ; do
	echo "$i" > ligne_en_cours

	scaffold=`cut -f1 ligne_en_cours` ; start=`cut -f2 ligne_en_cours` ; end=`cut -f3 ligne_en_cours` ;  gene=`cut -f4 ligne_en_cours` ; blast_start=`cut -f5 ligne_en_cours` ; blast_end=`cut -f6 ligne_en_cours`
	scaffold_length=`grep "^$scaffold	" $genome.fai | cut -f2`
	contig_start=30 ; contig_end=$((scaffold_length-30))



	samtools faidx $genome $scaffold:$blast_start-$blast_end > Region_Genomique_2.fa
	makeblastdb -in Region_Genomique_2.fa -dbtype nucl  
	samtools faidx ALL_functional_OR.prot $gene > Gene.fa



	tblastn -query Gene.fa -db Region_Genomique_2.fa -evalue 1e-20 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe sseq' > tblast_result

	cut -f12 tblast_result > Frame_detection
	cut -f13 tblast_result > Stop_detection

	number_frames=0
	Number_hsp=`cat Frame_detection | wc -l` 
	for i in `seq 2 $Number_hsp` ; do 
		j=$((i - 1))

		frame_first=`head -n $j Frame_detection| tail -1`
		frame_second=`head -n $i Frame_detection| tail -1`

		if [ $((frame_second - frame_first)) -ne '0' ] ; then
			number_frames=$((number_frames + 1))
		fi

	done

	frame=`cut -f12 Frame_detection | head -1`



	number_stops=`grep -o "\*" Stop_detection | wc -l`


	check_1=`if [ $end -gt $contig_end ] ; then echo "CONTIG_END" ; else echo "NA" ; fi`
	check_2=`if [ $start -lt $contig_start ] ; then echo "CONTIG_START" ; else echo "NA" ; fi`


	echo "$scaffold	$start	$end	$check_1	$check_2	$number_frames	$number_stops	$frame" >> Truncated_and_pseudo_coordinates.tsv

done



grep -v "CONTIG" Truncated_and_pseudo_coordinates.tsv | awk '($6 > 0 || $7 > 0)' > Pseudo_coordinates.tsv
grep -v "CONTIG" Truncated_and_pseudo_coordinates.tsv | awk '($6 == 0 && $7 == 0)' > Truncated_coordinates.tsv

grep "CONTIG" Truncated_and_pseudo_coordinates.tsv | awk '($6 > 0 || $7 > 0)' >> Pseudo_coordinates.tsv
grep "CONTIG" Truncated_and_pseudo_coordinates.tsv | awk '($6 == 0 && $7 == 0)' > Edges_coordinates.tsv



#Extract CDS sequences of truncated and pseudogenes and create corresponding proteins sequences file

for line in `cat Edges_coordinates.tsv` ; do echo "$line" > ligne_en_cours ; scaffold=`cut -f1 ligne_en_cours` ; start=`cut -f2 ligne_en_cours` ; end=`cut -f3 ligne_en_cours` ; if [ $frame -lt 0 ] ; then samtools faidx $genome $scaffold:$start-$end > temp.fa ; sed -i 's/:/-/g' temp.fa ; revseq temp.fa temp_rev.fa ; cat temp_rev.fa >> Edges_multifasta_temp.fa ; else samtools faidx $genome $scaffold:$start-$end >> Edges_multifasta_temp.fa ; fi ; done 
sed -i 's/Reversed://g' Edges_multifasta_temp.fa
sed -i 's/-/:/g' Edges_multifasta_temp.fa
cat Edges_multifasta_temp.fa Edge_multi_exon.fa > All_Edge_genes.fa
mv All_Edge_genes.fa Edges_multifasta_temp.fa 
sed -i 's/>/>E_/g' Edges_multifasta_temp.fa 

for line in `cat Pseudo_coordinates.tsv` ; do echo "$line" > ligne_en_cours ; scaffold=`cut -f1 ligne_en_cours` ; start=`cut -f2 ligne_en_cours` ; end=`cut -f3 ligne_en_cours` ; frame=`cut -f8 ligne_en_cours` ; if [ $frame -lt 0 ] ; then samtools faidx $genome $scaffold:$start-$end > temp.fa ; sed -i 's/:/-/g' temp.fa ; revseq temp.fa temp_rev.fa ; cat temp_rev.fa >> Pseudogenes_multifasta_temp.fa ; else samtools faidx $genome $scaffold:$start-$end >> Pseudogenes_multifasta_temp.fa ; fi ; done
sed -i 's/Reversed://g' Pseudogenes_multifasta_temp.fa
sed -i 's/-/:/g' Pseudogenes_multifasta_temp.fa
cat Pseudo_multi_exon.fa Pseudogenes_multifasta_temp.fa > All_Pseudo.fa
mv All_Pseudo.fa Pseudogenes_multifasta_temp.fa
sed -i 's/>/>P_/g' Pseudogenes_multifasta_temp.fa 

for line in `cat Truncated_coordinates.tsv` ; do echo "$line" > ligne_en_cours ; scaffold=`cut -f1 ligne_en_cours` ; start=`cut -f2 ligne_en_cours` ; end=`cut -f3 ligne_en_cours` ; frame=`cut -f8 ligne_en_cours` ; if [ $frame -lt 0 ] ; then samtools faidx $genome $scaffold:$start-$end > temp.fa ; sed -i 's/:/-/g' temp.fa ; revseq temp.fa temp_rev.fa ; cat temp_rev.fa >> Truncated_multifasta.fa ; else samtools faidx $genome $scaffold:$start-$end >> Truncated_multifasta.fa ; fi ; done
sed -i 's/Reversed://g' Truncated_multifasta.fa 
sed -i 's/-/:/g' Truncated_multifasta.fa 
cat Truncated_multi_exon.fa Truncated_multifasta.fa  > All_Truncated.fa
mv All_Truncated.fa Truncated_multifasta.fa 
sed -i 's/>/>T_/g' Truncated_multifasta.fa 



#Remove functional genes without 7tm (check with interproscan) and add it to truncated genes

IFS=$'\n'

rm -r temp/
rm temp

/home/policarpo/interproscan-5.35-74.0/interproscan.sh -i clear_Functionnal_ORs_multifasta_cdhit.prot -o 7tm_checking -f TSV

grep "IPR017452" 7tm_checking | cut -f1 | sort | uniq > Genes_without_gaps_in_7tm

grep ">" clear_Functionnal_ORs_multifasta_cdhit.prot | sed 's/>//g' > list_initiale

for gene in `cat list_initiale` ; do if grep -q "$gene" Genes_without_gaps_in_7tm ; then echo "clear" ; else echo "$gene" >> Genes_with_gaps_in_7tm ; fi ; done


for i in `cat Genes_without_gaps_in_7tm` ; do samtools faidx clear_Functionnal_ORs_multifasta_cdhit.fa $i ; done > RESULTS_Functionnal_genes.fa
for i in `cat Genes_with_gaps_in_7t` ; do samtools faidx clear_Functionnal_ORs_multifasta_cdhit.fa $i ; done >> Truncated_multifasta.fa



# DOnt forget to add multiexon functional genes to the final functional fasta file (and remove redundant sequences before)

/home/policarpo/cdhit-master/cd-hit -i Multiple_exon_functional_genes.fa -o Multiple_exon_functional_genes_cdhit.fa -c 1

cat Multiple_exon_functional_genes_cdhit.fa RESULTS_Functionnal_genes.fa > Temp_fasta
mv Temp_fasta RESULTS_Functionnal_genes.fa


#Counter the number of functional, edge, pseudogene, truncated and ambigous sequences

/home/policarpo/cdhit-master/cd-hit -i Edges_multifasta_temp.fa -o Edges_multifasta_cdhit.fa -c 1
/home/policarpo/cdhit-master/cd-hit -i Pseudogenes_multifasta_temp.fa -o Pseudogenes_multifasta_temp_cdhit.fa -c 1
/home/policarpo/cdhit-master/cd-hit -i Truncated_multifasta.fa -o Truncated_multifasta_cdhit.fa -c 1
/home/policarpo/cdhit-master/cd-hit -i Ambigous_sequences_multifasta.fa -o Ambigous_sequences_multifasta_cdhit.fa -c 1



x1=`grep -c ">" RESULTS_Functionnal_genes.fa`
x2=`grep -c ">" Edges_multifasta_cdhit.fa` 
x3=`grep -c ">" Pseudogenes_multifasta_temp_cdhit.fa`
x4=`grep -c ">" Truncated_multifasta_cdhit.fa` 
x5=`grep -c ">" Ambigous_sequences_multifasta_cdhit.fa` 


echo "$x1	$x2	$x3	$x4	$x5" > Count_ORs.tsv











