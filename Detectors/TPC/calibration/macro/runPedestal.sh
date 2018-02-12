#!/bin/bash

if [ $# -lt 3 ]; then
  echo "usage: runPedestal <fileInfo> <pedestalFile> <nevents>"
fi

fileInfo=$1
pedestalFile=$2
nevents=$3

adcMin=0
adcMax=1100
slowFit=0

if [ $# -ge 5 ]; then
  adcMin=$4
  adcMax=$5
  slowFit=$6
fi

cmd="root.exe -b -q -l -n -x $O2_SRC/Detectors/TPC/calibration/macro/runPedestal.C'(\"$fileInfo\",\"$pedestalFile\", $nevents, $adcMin, $adcMax, $slowFit)'"
echo "running: $cmd"
eval $cmd
