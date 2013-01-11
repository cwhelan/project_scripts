#!/bin/bash

set -e
set -u

LABEL=$1
RANKED_GENE_LIST=$2
GENESET=$3

cat <<EOF
java -Xmx4g -cp /u0/dbase/cw/software/lib/gsea2-2.08.jar \
    xtools.gsea.GseaPreranked \
    -gmx $GENESET \
    -collapse false -mode Max_probe -norm meandiv -nperm 1000 \
    -rnk ${RANKED_GENE_LIST} \
    -scoring_scheme weighted -rpt_label $LABEL -include_only_symbols true \
    -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report true \
    -out ${RANKED_GENE_LIST}.gsea.out -gui false
EOF

java -Xmx4g -cp /u0/dbase/cw/software/lib/gsea2-2.08.jar \
    xtools.gsea.GseaPreranked \
    -gmx $GENESET \
    -collapse false -mode Max_probe -norm meandiv -nperm 1000 \
    -rnk ${RANKED_GENE_LIST} \
    -scoring_scheme weighted -rpt_label $LABEL -include_only_symbols true \
    -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report true \
    -out ${RANKED_GENE_LIST}.gsea.out -gui false