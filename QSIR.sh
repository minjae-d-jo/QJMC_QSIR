#!/bin/bash

endTime=$1
N=$2
omega=$3
nRepeat=$4
uid=$5

[[ $# -ne 5 ]] && echo "Error: You need to give 5 arguments." && exit 1

params=${endTime}_${N}_${omega}_${nRepeat}_${uid}
logFile=Log/QJMC_$params.table
outputFile=Result/QJMC_$params.table
#outputFile2=Result/DP_Order_1D_$params.table

exec 2>> $logFile

echo "$(date +%Y-%m-%d_%H:%M:%S): start $$ $(hostname)" >> $logFile

Bin/QSIR_${N} $endTime $omega $nRepeat >> $outputFile #2>> $outputFile2
ret=$?;
echo -n "$(date +%Y-%m-%d_%H:%M:%S): " >> $logFile
if [ $ret -eq 0 ]; then
	echo "normal exit"
else
	echo "abnormal exit: $ret"
fi >> $logFile
