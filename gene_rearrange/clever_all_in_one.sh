#!/bin/bash

set -e
set -u

BAM_FILE=$1
REF=$2

WORKING_DIR=`dirname $BAM_FILE`
FILE_NAME=`basename $BAM_FILE`
FILE_PREFIX=`echo $FILE_NAME | awk -F. '{ print $1 }'`

OUTPUT_DIR=${FILE_PREFIX}_clever_all_in_one

echo "clever-all-in-one -B -T 4 -I $BAM_FILE $REF $OUTPUT_DIR"
clever-all-in-one -k -B -T 4 -I  $BAM_FILE $REF $OUTPUT_DIR
