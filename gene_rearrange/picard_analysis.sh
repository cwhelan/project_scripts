#!/bin/bash

echo "Validate Sam..."
java -Xmx2g -jar /g/whelanch/software/picard-tools-1.72/ValidateSamFile.jar I=$1
echo "Collect alignment summary metrics.."
java -Xmx2g -jar /g/whelanch/software/picard-tools-1.72/CollectAlignmentSummaryMetrics.jar I=$1 O=$2
echo "Collect insert size metrics.."
java -Xmx2g -jar /g/whelanch/software/picard-tools-1.72/CollectInsertSizeMetrics.jar I=$1 O=$3 HISTOGRAM_FILE=$4
