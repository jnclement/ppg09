#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.521
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
export TESTINSTALL=/sphenix/user/jocl/projects/testinstall
echo $LD_LIBRARY_PATH
echo $PATH
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
N=$(( 100 * $1 ))
tail -n +$N seglist_data.list | head -n 100 > thelist.list
mkdir input
for i in {0..99}; do
    n=$(( $i + 1 ))
    FILENAME=`sed -n "$n"p thelist.list`
    cp $FILENAME input/inputfile_$(( $N + $i )).root
done

mkdir output

root -l -b -q "analyze_segment_data.C($N,100)"
/sphenix/user/jocl/projects/ppg09/ana_output/dat -p
cp output/* /sphenix/user/jocl/projects/ppg09/ana_output/dat
