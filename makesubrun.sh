#!/bin/bash

if [ $# -lt 2 ]; then
    echo "need TAG and ISDAT args (tag should be data if data, and jetxx otherwise)"
    exit 1
fi
TAG=$1
ISDAT=$2
TOSUB=$(( `cat seglist_$TAG.list | wc -l` + 99))
TOSUB=$(( $TOSUB / 100 ))
BASENAME="condor_${TAG}"
PREFIX="."
SUBNAME="${BASENAME}.sub"

EXE="sim_condor.sh"
if [[ $ISDAT -ne 0 ]]; then
    EXE="data_condor.sh"
fi
echo "executable = ${EXE}" > $PREFIX/$SUBNAME
echo "arguments = \$(Process) ${TAG}" >> $PREFIX/$SUBNAME
echo "output = output/out/output_${BASENAME}_\$(Process).out" >> $PREFIX/$SUBNAME
echo "request_memory                = 4GB" >> $PREFIX/$SUBNAME
echo "error = output/out/output_${BASENAME}_\$(Process).out" >> $PREFIX/$SUBNAME
echo "log = /tmp/jocl_${BASENAME}.log" >> $PREFIX/$SUBNAME
echo "queue ${TOSUB}" >> $PREFIX/$SUBNAME

condor_submit $PREFIX/$SUBNAME
