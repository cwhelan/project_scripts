#!/bin/bash

set -e
set -u

BP_FILE=$1
BP_NAME=$2
FEATURE_FILE=$3
FEATURE_NAME=$4

DIR_NAME=`echo "${BP_NAME}_${FEATURE_NAME}" | sed "s/ /_/g"`
echo "Setting up dir $DIR_NAME"

ESCAPED_BP_FILE=`echo $BP_FILE | sed "s/\//\\\//g"`
echo "Escaped BP File: $ESCAPED_BP_FILE"
ESCAPED_FEATURE_FILE=`echo $FEATURE_FILE | sed "s/\//\\\//g"`
echo "Escaped Feature File: $ESCAPED_FEATURE_FILE"

mkdir $DIR_NAME
cat permutation_tests.desc.tmpl | sed "s:BP_FILE:${BP_FILE}:" | sed "s/BP_NAME/${BP_NAME}/" | sed "s:FEATURE_FILE:${FEATURE_FILE}:" | sed "s/FEATURE_NAME/${FEATURE_NAME}/" | sed "s/DIR_NAME/${DIR_NAME}/" > $DIR_NAME/permutation_tests.desc


