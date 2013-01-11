#!/bin/bash

GENE_SET=$1

#/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh normalGenes_${GENE_SET} normalGene.rnk ${GENE_SET} &> normalGenes_${GENE_SET}.log &
/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh cancerGenes_${GENE_SET} cancerGene.rnk ${GENE_SET} &> cancerGenes_${GENE_SET}.log &
# /u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh intermediateGenes_${GENE_SET} intermediateGene.rnk ${GENE_SET} &> intermediateGenes_${GENE_SET}.log &
#/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh normalCancerGeneDiff_${GENE_SET} normalCancerDiff.rnk ${GENE_SET} &> normalCancerGeneDiff_${GENE_SET}.log &
#/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh normalIntermediateGeneDiff_${GENE_SET} normalIntermediateDiff.rnk ${GENE_SET} &> normalIntermediateGeneDiff_${GENE_SET}.log &
#/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh intermediateCancerGeneDiff_${GENE_SET} intermediateCancerDiff.rnk ${GENE_SET} &> intermediateCancerGeneDiff_${GENE_SET}.log &
#/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh normalPromoters_${GENE_SET} normalPromoter.rnk ${GENE_SET} &> normalPromoter_${GENE_SET}.log &
/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh cancerPromoters_${GENE_SET} cancerPromoter.rnk ${GENE_SET} &> cancerPromoter_${GENE_SET}.log &
/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh intermediatePromoters_${GENE_SET} intermediatePromoter.rnk ${GENE_SET} &> intermediatePromoter_${GENE_SET}.log &
#/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh normalCancerPromoters_${GENE_SET} normalCancerPromoterDiff.rnk ${GENE_SET} &> normalCancerPromoter_${GENE_SET}.log &
/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh intermediateCancerPromoters_${GENE_SET} intermediateCancerPromoterDiff.rnk ${GENE_SET} &> intermediateCancerPromoter_${GENE_SET}.log &
#/u0/dbase/cw/project_scripts/methyl-analysis/gsea.sh normalIntermediatePromoters_${GENE_SET} normalIntermediatePromoterDiff.rnk ${GENE_SET} &> normalIntermediatePromoter_${GENE_SET}.log &