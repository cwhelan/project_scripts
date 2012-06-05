#!/bin/bash

echo "mrfast --search $1 --seq $2/$4.$3 -o $2/$3_se.out"
mrfast --search $1 --seq $2/$4.$3 -o $2/$3_se.out
