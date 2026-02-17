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
cp /sphenix/user/jocl/projects/ppg09/seglist_$2.list .
N=$(( 20 * $1 ))
N=$(( $N + 1 ))
tail -n +$N seglist_$2.list | head -n 20 > thelist.list
cat thelist.list
mkdir input
cp /sphenix/user/jocl/projects/ppg09/input_zvertexreweight.root  .
cp /sphenix/user/jocl/projects/ppg09/output_jetefficiency.root  .
cp /sphenix/user/jocl/projects/ppg09/output_mbdefficiency.root .
cp /sphenix/user/jocl/projects/ppg09/unfold_Def.h .
cp /sphenix/user/jocl/projects/ppg09/analyze_segment_sim.C .
cp /sphenix/user/jocl/projects/ppg09/output_reweightfunction_r04.root .
for i in {0..19}; do
    n=$(( $i + 1 ))
    FILENAME=`sed -n "$n"p thelist.list`
    cp $FILENAME input/inputfile_$i.root
done

mkdir output

root -l -b -q 'analyze_segment_sim.C("'$2'",'$N',20,'$3')'
mkdir -p /sphenix/user/jocl/projects/ppg09/ana_output/$2
cp output/* /sphenix/user/jocl/projects/ppg09/ana_output/$2
