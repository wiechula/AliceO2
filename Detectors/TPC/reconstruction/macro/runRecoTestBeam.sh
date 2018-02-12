#!/bin/bash

if [ $# -lt 3 ]; then
  echo "usage: runTracking <fileInfo> <pedestalFile> <nevents>"
fi

fileInfo=$1
pedestalFile=$2
nevents=$3

clusterFile="clusters.root"
trackFile="tracks.root"

script=$(readlink -f $0)
macroDir=$(dirname $script)
findClusters=${macroDir}/RawClusterFinder.C
runTracking=${macroDir}/runCATrackingTestBeam.C
addInclude=${macroDir}/addInclude.C
cherenkovFile=$PWD/cherenkov.txt

# ===| find raw clusters |======================================================
if [ ! -f ${clusterFile} ]; then
  cmd="root.exe -b -q -x -l  ${addInclude} $findClusters'+(\"$fileInfo\",\"$pedestalFile\",\"$clusterFile\",$nevents)'"
  echo $cmd
  eval $cmd
else
   echo "${clusterFile} already exists, delete it first to rerun."
     echo " You should then also delete ${trackFile}."
fi

# ===| run reconstruction |=====================================================
if [ ! -f ${trackFile} ]; then
  if [ -f ${clusterFile} ]; then
    cmd="root.exe -b -q -x -l  ${addInclude} $runTracking'+(\"$clusterFile\",\"$trackFile\")'"
    echo $cmd
    eval $cmd
  else
    echo "${clusterFile} does not exists!"
  fi
else
   echo "${trackFile} already exists, delete it first to rerun."
fi

