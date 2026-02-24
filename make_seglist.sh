#!/bin/bash

#rm seglist_data.list
#for rn in `cat listrunnumber.txt`; do
#    echo $rn
#    ls -v /sphenix/tg/tg01/jets/jocl/chi2/$rn/*20260216*chi2file* >> seglist_data.list #ls /sphenix/tg/tg01/jets/jocl/chi2/$rn/*20260116_dat*chi2file* >> seglist_data.list
#done

for stype in mb; do # jet5 jet12 jet20 jet30 jet40 jet50 jet60; do
    find /sphenix/tg/tg01/jets/jocl/chi2/ -type f -name "*20260219_${stype}_*chi2file*" -print | sort -V > seglist_$stype.list
done
