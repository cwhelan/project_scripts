#!/bin/bash

BISMARK_OUT=$1

echo "methylation_extractor -p --comprehensive --no_overlap --report $1"
methylation_extractor -p --comprehensive --no_overlap --report $1
