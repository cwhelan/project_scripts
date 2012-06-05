#!/bin/bash

echo "mrfast -b --search $1 --pe --seq1 $2/$4.$3 --seq2 $2/$5.$3 --discordant-vh --min $6 --max $7 -o $2/$3.out"
mrfast -b --search $1 --pe --seq1 $2/$4.$3 --seq2 $2/$5.$3 --discordant-vh --min $6 --max $7 -o $2/$3.out 
