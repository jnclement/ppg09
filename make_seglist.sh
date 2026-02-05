#!/bin/bash

#rm seglist_data.list
#for rn in `cat listrunnumber.txt`; do
#    echo $rn
#    ls /sphenix/tg/tg01/jets/jocl/chi2/$rn/*20250130*chi2file* >> seglist_data.list #ls /sphenix/tg/tg01/jets/jocl/chi2/$rn/*20260116_dat*chi2file* >> seglist_data.list
#done

for stype in jet5 jet10 jet15 jet20 jet30 jet50 jet70; do
    rm seglist_$stype.list
    for i in {0..499}; do
	echo $stype $i
	ls /sphenix/tg/tg01/jets/jocl/chi2/$i/*20260202_$stype*chi2file* >> seglist_$stype.list
    done
done

rm seglist_mb.list
for i in {0..1999}; do
    echo $i
    ls /sphenix/tg/tg01/jets/jocl/chi2/$i/*20260202_mb*chi2file* >> seglist_mb.list
done
