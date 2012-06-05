#!/bin/bash

for ((i=17; i<=450; i++))
do
	printf "processing split file number %3.3d\n" $i
	echo "qsub -cwd -V -l q=*kiwi* -o job_out/ -e job_out/ ./mrfast_align_file_pair_single.sh /l2/users/whelanch/1000genomes/ref/human_g1k_v37.fasta /l2/users/whelanch/gene_rearrange/clean_splits `printf "%3.3d" $i` s_3_1_rc s_3_2_rc 3000 8000"
	qsub -cwd -V -l q=*kiwi* -o job_out/ -e job_out/ ./mrfast_align_file_pair_single.sh /l2/users/whelanch/1000genomes/ref/human_g1k_v37.fasta /l2/users/whelanch/gene_rearrange/clean_splits `printf "%3.3d" $i` s_3_1_rc s_3_2_rc 3000 8000
	echo "sleeping..."
 	sleep 50m

done
