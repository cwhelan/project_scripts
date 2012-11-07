#!/bin/bash

set -e
set -u

PINDEL_D_FILE=$1

cat $PINDEL_D_FILE | grep ChrID | awk '{OFS="\t"; print $8,$10,$11,$25}'