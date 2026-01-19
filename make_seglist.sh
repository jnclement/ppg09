#!/bin/bash

rm seglist_data.list
for rn in `cat listrunnumber.txt`; do
    ls /sphenix/tg/tg01/jets/jocl/chi2/$rn/*20260116_dat*chi2file* >> seglist_data.list
done

exit 0

for stype in jet5 jet10 jet15 jet20 jet30 jet50 jet70; do
    rm seglist_$stype.list
    for i in {0..10000}; do
	ls /sphenix/tg/tg01/jetes/jocl/chi2/$i/*20260116_$stype*chi2file* >> seglist_$stype.list
    done
done

rm seglist_mb.list
for i in {0..200000}; do
    ls /sphenix/tg/tg01/jetes/jocl/chi2/$i/*20260116_mb*chi2file* >> seglist_mb.list
done
