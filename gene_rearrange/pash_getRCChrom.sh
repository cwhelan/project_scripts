#!/bin/bash

REF_FILE=`basename $1`
REF_DIR=`dirname $1`
getRCChrom.rb $1 /l2/users/whelanch/gene_rearrange/pash/$REF_FILE.dnameth.fa

