#!/bin/bash

if [ $# -lt 1 ]; then
    echo "need type argument (jetX or mb)"
    exit 1
fi

NPROC=100
if [ "$1" == "mb" ]; then
    NPROC=400
fi
NPROC=2
rm condor_ana_$1.sub
echo "executable = sim_condor.sh" >> condor_ana_$1.sub
echo "arguments = \$(Process) ${1}" >> condor_ana_$1.sub
echo "output = output/out/output_condor_${1}_\$(Process).out" >> condor_ana_$1.sub
echo "request_memory                = 4GB" >> condor_ana_$1.sub
echo "error = output/out/output_condor_${1}_\$(Process).out" >> condor_ana_$1.sub
echo "log = /tmp/jocl_condor_${1}.log" >> condor_ana_$1.sub
echo "queue ${NPROC}" >> condor_ana_$1.sub
condor_submit condor_ana_$1.sub
