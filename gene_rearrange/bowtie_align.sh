
echo "bowtie -q --solexa1.3-quals -n 3 -I 3000 -X 8000 --rf -S -p 8 --chunkmbs 1024 /l2/users/whelanch/1000genomes/ref/hgv37_bowtie -1 /l2/users/whelanch/gene_rearrange/s_3_1_sequence.txt -2 /l2/users/whelanch/gene_rearrange/s_3_2_sequence.txt /l2/users/whelanch/gene_rearrange/bowtie_rf_3.out"
bowtie -q --solexa1.3-quals -n 3 -I 3000 -X 8000 --rf -S -p 8 --chunkmbs 1024 /l2/users/whelanch/1000genomes/ref/hgv37_bowtie -1 /l2/users/whelanch/gene_rearrange/s_3_1_sequence.txt -2 /l2/users/whelanch/gene_rearrange/s_3_2_sequence.txt /l2/users/whelanch/gene_rearrange/bowtie_rf_3.out
