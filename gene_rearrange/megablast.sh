#!/bin/bash

for f in tier1*
do
	/g/whelanch/software/bin/megablast -d hg19 -i $f -W 8 -G 6 -E 4 -F F -q -3 -r 2 -D 3 > $f.mb
done
