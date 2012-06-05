#!/bin/bash

f1=`ls s_3_1*`
f2=`ls s_3_2*`
echo "mrfast --search $1 --pe --seq1 $f1 --seq2 $f2 --discordant-vh --min $2 --max $3 -o mrfast.out"

mrfast --search $1 --pe --seq1 $f1 --seq2 $f2 --discordant-vh --min $2 --max $3 -o mrfast.out
