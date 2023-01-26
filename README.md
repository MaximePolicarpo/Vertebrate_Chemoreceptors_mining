# Vertebrate_Chemoreceptors_mining

The main pipelines scripts are called OR_Finder.sh, TAAR_Finder.sh, V1R_Finder.sh, V2R_Finder.sh, T1R_Finder.sh, T2R_Finder.sh

1) Before running these pipelines, create a conda environment called "olfactory" and install seqkit, scikit-learn and tmhmm : 

conda create -n olfactory
conda activate olfactory
conda install -c "bioconda/label/cf201901" seqkit
conda install -c anaconda scikit-learn
conda install -c dansondergaard tmhmm.py


2) If you are not under a slurm environment, then you should also have the following programs installed on your machine : 

R v4.2.0
BLAST v2.12.0
EMBOSS v6.2.0
SAMtools v1.15
MAFFT v7.467
IQ-TREE v2.0
Python v3.9.5
FASTX-Toolkit v0.0.14


3a) If you are under a slurm environment, make sure that you have a qos named 6hours. Otherwise, replace the qos name in the lines beginning with "sbatch" in the main .sh scripts [OR_Finder.sh, TAAR_Finder.sh, V1R_Finder.sh, V2R_Finder.sh, T1R_Finder.sh, T2R_Finder.sh] (for example line 222 in TAAR_Finder.sh). 6hours is the optimal running time for these sbatch commands. 


3b) Again, if you are not under a slurm environment, then replace lines beginning by sbatch with nohup

For example in TAAR_Finder.sh : 
Replace

sbatch -W -c 4 --qos=6hours --wrap="$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --ryo '%tcs' --minintron 50 --maxintron $maximum_intron_length Splitted_db/$file_name TAAR_best_hits_regions.fa > Exonerate_raw_results_folder/$file_name.exo.rslt ; sleep 10" &

by 

nohup $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --ryo '%tcs' --minintron 50 --maxintron $maximum_intron_length Splitted_db/$file_name TAAR_best_hits_regions.fa > Exonerate_raw_results_folder/$file_name.exo.rslt &


4) OR genes

If you want to extract OR genes from a genome : 

sbatch OR_Finder.sh $genome_file_name $OR_database $GPCR_database $Script_folders_location $intron_sizes $thread_number $phobius_tmhmm $max_nb_exons

- $genome_file_name : You genome fasta file
- $OR_database : Database containing protein sequences known OR genes. You can either create one of your choice or use the one provided here (Database_OR/Database_Vertebrate_OR_cdhit_80.prot)
- $GPCR_database : Blast database of non-chemoreceptors GPCR proteins + chemoreceptors proteins. One should use the database provided (GPCR_plus_Chemoreceptors_vertebrates.prot)
- $Script_folders_location : Full path the the folder containing accessory OR_finder scripts (Database_OR/Scripts_2022/)
- $intron_sizes : Maximum intron size for multi-exon TAAR genes. (optimal : 25000)
- $thread_number : Number of thread that will be used by IQ-TREE and BLAST
- $phobius_tmhmm : Set to TRUE or FALSE. If TRUE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed into the same file as pseudogenes and truncated genes. If FALSE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed in the result file together with genes with a predicted 7tm domain. 
- $max_nb_exons : Maximum number of exons a OR gene should have. The opitmal is to set 2 for tetrapod OR genes, and 4 for ray-finned fishes.



5) TAAR genes

If you want to extract TAAR genes from a genome : 


sbatch TAAR_Finder.sh $genome_file_name $TAAR_database $GPCR_database $Script_folders_location $intron_sizes $thread_number $Exon_mode $evalue $phobius_tmhmm

- $genome_file_name : You genome fasta file
- $TAAR_database : Database containing protein sequences known OR genes. You can either create one of your choice or use the one provided here (Database_TAAR/TAAR_plus_TAARL_database_reformat_cdhit_80.prot)
- $GPCR_database : Blast database of non-chemoreceptors GPCR proteins + chemoreceptors proteins. One should use the database provided (GPCR_plus_Chemoreceptors_vertebrates.prot)
- $Script_folders_location : Full path the the folder containing accessory TAAR_finder scripts (Database_TAAR/Scripts_2022/)
- $intron_sizes : Maximum intron size for multi-exon TAAR genes. (optimal : 25000)
- $thread_number : Number of thread that will be used by IQ-TREE and BLAST
- $Exon_mode : if TRUE, then the pipeline will look for the presence of multiple exons TAAR genes. If FALSE, then only one exon genes will be retrieved but the pipeline will run much faster
- $evalue : evalue for the initial tblastn of known TAAR genes against the genome. Optimal : 1e-5
- $phobius_tmhmm : Set to TRUE or FALSE. If TRUE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed into the same file as pseudogenes and truncated genes. If FALSE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed in the result file together with genes with a predicted 7tm domain.


5) V1R genes 

If you want to extract V1R genes from a genome : 

sbatch V1R_Finder.sh $genome_file_name $V1R_database $GPCR_database $Script_folders_location $intron_sizes $thread_number $Exon_mode $evalue $phobius_tmhmm


- $genome_file_name : You genome fasta file
- $V1R_database : Database containing protein sequences known V1R genes. You can either create one of your choice or use the one provided here (Database_V1R/Database_2022_V1R_vertebrates_cdhit_80.prot)
- $GPCR_database : Blast database of non-chemoreceptors GPCR proteins + chemoreceptors proteins. One should use the database provided (GPCR_plus_Chemoreceptors_vertebrates.prot)
- $Script_folders_location : Full path the the folder containing accessory V1R_finder scripts (Database_V1R/Scripts_2022/)
- $intron_sizes : Maximum intron size for multi-exon V1R genes. (optimal : 25000)
- $thread_number : Number of thread that will be used by IQ-TREE and BLAST
- $Exon_mode : if TRUE, then the pipeline will look for the presence of multiple exons V1R genes. If FALSE, then only one exon genes will be retrieved but the pipeline will run much faster. One MUST set this to TRUE for ray-finned genomes. 
- $evalue : evalue for the initial tblastn of known TAAR genes against the genome. Optimal : 1e-5
- $phobius_tmhmm : Set to TRUE or FALSE. If TRUE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed into the same file as pseudogenes and truncated genes. If FALSE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed in the result file together with genes with a predicted 7tm domain.



6) V2R genes 

If you want to extract V2R genes from a genome : 


sbatch V2R_Finder.sh $genome_file_name $V2R_database $GPCR_database $Script_folders_location $intron_sizes $thread_number $phobius_tmhmm

- $genome_file_name : You genome fasta file
- $V2R_database : Database containing protein sequences known V2R genes. You can either create one of your choice or use the one provided here (Database_V2R/Database_2022_V2R_vertebrates_cdhit_80.prot)
- $GPCR_database : Blast database of non-chemoreceptors GPCR proteins + chemoreceptors proteins. One should use the database provided (GPCR_plus_Chemoreceptors_vertebrates.prot)
- $Script_folders_location : Full path the the folder containing accessory V2R_finder scripts (Database_V2R/Scripts_2022/)
- $intron_sizes : Maximum intron size for V2R genes. (Set between 15,000 ad 40,000 depending on the species)
- $thread_number : Number of thread that will be used by IQ-TREE and BLAST
- $phobius_tmhmm : Set to TRUE or FALSE. If TRUE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed into the same file as pseudogenes and truncated genes. If FALSE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed in the result file together with genes with a predicted 7tm domain.




7) T1R genes 

If you want to extract T1R genes from a genome : 


sbatch T1R_Finder.sh $genome_file_name $T1R_database $GPCR_database $Script_folders_location $intron_sizes $thread_number $phobius_tmhmm

- $genome_file_name : You genome fasta file
- $T1R_database : Database containing protein sequences known T1R genes. You can either create one of your choice or use the one provided here (Database_T1R/Database_2022_T1R_vertebrates_cdhit_80.prot)
- $GPCR_database : Blast database of non-chemoreceptors GPCR proteins + chemoreceptors proteins. One should use the database provided (GPCR_plus_Chemoreceptors_vertebrates.prot)
- $Script_folders_location : Full path the the folder containing accessory T1R_finder scripts (Database_T1R/Scripts_2022/)
- $intron_sizes : Maximum intron size for T1R genes. (Set between 15,000 ad 40,000 depending on the species)
- $thread_number : Number of thread that will be used by IQ-TREE and BLAST
- $phobius_tmhmm : Set to TRUE or FALSE. If TRUE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed into the same file as pseudogenes and truncated genes. If FALSE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed in the result file together with genes with a predicted 7tm domain.


8) T2R genes


sbatch T2R_Finder.sh $genome_file_name $T2R_database $GPCR_database $Script_folders_location $intron_sizes $thread_number $evalue $phobius_tmhmm


 - $genome_file_name : You genome fasta file
- $T2R_database : Database containing protein sequences known T2R genes. You can either create one of your choice or use the one provided here (Database_T2R/Vertebrates_T2R_db_cdhit_70.prot)
- $GPCR_database : Blast database of non-chemoreceptors GPCR proteins + chemoreceptors proteins. One should use the database provided (GPCR_plus_Chemoreceptors_vertebrates.prot)
- $Script_folders_location : Full path the the folder containing accessory T2R_finder scripts (Database_T2R/Scripts_2022/)
- $intron_sizes : Maximum intron size for T2R genes (Deprecated, you can put a random number. I never found any multiple-exon T2R genes)
- $thread_number : Number of thread that will be used by IQ-TREE and BLAST
- $evalue : evalue for the initial tblastn of known T2R genes against the genome. Optimal : 1e-5
- $phobius_tmhmm : Set to TRUE or FALSE. If TRUE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed into the same file as pseudogenes and truncated genes. If FALSE, then genes without a 7tm domain predicted by phobius and/or tmhmm will be placed in the result file together with genes with a predicted 7tm domain.



9) Example

Lets say I want to extract T2R genes from the Danio rerio genome. First I download the zebrafish genome and unzip it (fasta):

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/035/GCA_000002035.4_GRCz11/GCA_000002035.4_GRCz11_genomic.fna.gz
gzip -d GCA_000002035.4_GRCz11_genomic.fna.gz


Then run the T2R pipeline, with a tblastn evalue of 1e-5 and 20 threads : 

sbatch T2R_Finder.sh GCA_000002035.4_GRCz11_genomic.fna ./Database_T2R/Vertebrates_T2R_db_cdhit_70.prot GPCR_plus_Chemoreceptors_vertebrates.prot ./Database_T2R/Scripts_2022/ 25000 20 1e-5 TRUE


