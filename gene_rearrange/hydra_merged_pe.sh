#!/bin/bash

set -e
set -u

BAM_FILE=$1
INSERT_SIZE_STAT_FILE=$2

WORKING_DIR=`dirname $BAM_FILE`
FILE_NAME=`basename $BAM_FILE`
FILE_PREFIX=`echo $FILE_NAME | awk -F. '{ print $1 }'`

MEDIAN_ISIZE=`cat $INSERT_SIZE_STAT_FILE | grep SI | awk '{print $3}'`
MAD_ISIZE=`cat $INSERT_SIZE_STAT_FILE | grep SI | awk '{print $6}'`


MAX_CONCORDANT=`echo "($MEDIAN_ISIZE + 4 * $MAD_ISIZE + 0.5) / 1" | bc`
MIN_CONCORDANT=`echo "($MEDIAN_ISIZE - 4 * $MAD_ISIZE + 0.5) / 1" | bc`

echo "bamToBed -i $BAM_FILE -tag NM pairDiscordants.py -i stdin -m hydra -y $MIN_CONCORDANT -z $MAX_CONCORDANT > $WORKING_DIR/$FILE_PREFIX.bedpe"
bamToBed -ed -i $BAM_FILE | pairDiscordants.py -i stdin -m hydra -y $MIN_CONCORDANT -z $MAX_CONCORDANT > $WORKING_DIR/$FILE_PREFIX.bedpe

#echo "bamToBed -i $BAM_FILE -tag NM | pairDiscordants.py -i stdin -m hydra -y 2574 -z 7278 > $WORKING_DIR/$FILE_PREFIX.bedpe"
#bamToBed -i $BAM_FILE -tag NM | \
#    pairDiscordants.py -i stdin -m hydra -y 2574 -z 7278 > $WORKING_DIR/$FILE_PREFIX.bedpe

echo "dedupDiscordants.py -i $WORKING_DIR/$FILE_PREFIX.bedpe -s 3 > $WORKING_DIR/$FILE_PREFIX.dedup.bedpe"
dedupDiscordants.py -i $WORKING_DIR/$FILE_PREFIX.bedpe -s 3 > $WORKING_DIR/$FILE_PREFIX.dedup.bedpe

MAX_ALLOWABLE_LEN_DIF=`echo "(10 * $MAD_ISIZE + 0.5) / 1" | bc`
MAX_ALLOWABLE_NON_OVERLAP=`echo "(10 * $MEDIAN_ISIZE + 20 * $MAD_ISIZE + 0.5) / 1" | bc`

echo "hydra -in $WORKING_DIR/$FILE_PREFIX.dedup.bedpe -out $WORKING_DIR/${FILE_PREFIX}_all.breaks -mld $MAX_ALLOWABLE_LEN_DIF -mno $MAX_ALLOWABLE_NON_OVERLAP -li -use all"
hydra -in $WORKING_DIR/$FILE_PREFIX.dedup.bedpe -out $WORKING_DIR/${FILE_PREFIX}_all.breaks -mld $MAX_ALLOWABLE_LEN_DIF -mno $MAX_ALLOWABLE_NON_OVERLAP -li -use all
echo "hydra -in $WORKING_DIR/$FILE_PREFIX.dedup.bedpe -out $WORKING_DIR/${FILE_PREFIX}_ms3_all.breaks -mld $MAX_ALLOWABLE_LEN_DIF -mno $MAX_ALLOWABLE_NON_OVERLAP -ms 3 -li -use all"
hydra -in $WORKING_DIR/$FILE_PREFIX.dedup.bedpe -out $WORKING_DIR/${FILE_PREFIX}_ms3_all.breaks -mld $MAX_ALLOWABLE_LEN_DIF -mno $MAX_ALLOWABLE_NON_OVERLAP -ms 3 -li -use all
