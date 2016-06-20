#! /bin/bash
MainFile=$1
##Loop over samples
cd $MainFile

for ((i=1; i<=10  ; i++)); do
cd S$i
    python ~/pizza-9Oct15/src/pizza.py -f ../../V2CountOrientationCorr.py  
    wait

cd ..
done

