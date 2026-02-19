#!/bin/bash

if [ $# -lt 3 ]; then
    echo "need TAG, ISDAT, RADIUS args (tag should be data if data, and jetxx otherwise)"
    exit 1
fi
TAG=$1
ISDAT=$2
RAD=$3
TOSUB=$(( `cat seglist_$TAG.list | wc -l` + 19))
TOSUB=$(( $TOSUB / 20 ))
BASENAME="condor_${TAG}"
PREFIX="."
SUBNAME="${BASENAME}.sub"

EXE="sim_condor.sh"
if [[ $ISDAT -ne 0 ]]; then
    EXE="data_condor.sh"
    TOSUB=`cat listrunnumber.txt | wc -l`
fi
#TOSUB=1
echo "executable = ${EXE}" > $PREFIX/$SUBNAME
echo "arguments = \$(Process) ${TAG} ${RAD}" >> $PREFIX/$SUBNAME
echo "output = output/out/output_${BASENAME}_\$(Process).out" >> $PREFIX/$SUBNAME
echo "request_memory                = 4GB" >> $PREFIX/$SUBNAME
echo "error = output/out/output_${BASENAME}_\$(Process).out" >> $PREFIX/$SUBNAME
echo "log = /tmp/jocl_${BASENAME}.log" >> $PREFIX/$SUBNAME
echo "queue ${TOSUB}" >> $PREFIX/$SUBNAME

condor_submit $PREFIX/$SUBNAME
