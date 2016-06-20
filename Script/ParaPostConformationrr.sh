MainFile="$1"
##Loop over samples
cd $MainFile

for ((i=1; i<=10  ; i++)); do
cd S$i
    python ~/0_SourceCode_WithVersionControl/pizza-2Jul14/src/pizza.py -f ../../Countrr.py  
    wait

cd ..
done

