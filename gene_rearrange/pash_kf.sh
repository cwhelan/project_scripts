#!/bin/bash

REF_FILE=`basename $1`
REF_DIR=`dirname $1`
keyFreq.exe -o /l2/users/whelanch/gene_rearrange/pash/$REF_FILE.12.18.kf -p 111010110100110111 $1
keyFreq.exe -h -o /l2/users/whelanch/gene_rearrange/pash/$REF_FILE.12.18.kf.h -p 111010110100110111 $1


